function [Tx,Fx,Ty,Fy,Th,Fh,Xmat,tauxInt,tauyInt,etaInt,HeInt]=...
    uvhElementKernel(MUA,conn,...
    bv,hv,uv,vv,AGlenv,nv,Cv,mv,qv,mukv,V0v,h0v,u0v,v0v,as0v,ab0v,as1v,ab1v,dadhv,Bv,Sv,rhov,...
    LSFMaskv,hBCsMaskv,...
    uov,vov,Cov,mov,uav,vav,Cav,mav,...
    CtrlVar,rhow,g,Ronly,ca,sa,dt)

%% Element-based uvh assembly kernel (25 Oct 2026), compiled with MATLAB Coder, see README.md in this folder
%
% CODEGEN VERSION: options and parameters in CtrlVar are run-time inputs; only Ronly and the derivative flags are constants.
% Inputs are the nodal vectors (Nnodes x 1) and the connectivity; the nodal values of each element are gathered here. The element-matrix
% values are returned element-major in Xmat (9*nod^2 x Nele) and summed into the nonzeros by GatherSumCSC (see uvhMatrixAssemblyCompiled.m). CtrlVar.BasalDrag.CalculateDerivatives and
% CtrlVar.EffectiveViscosity.CalculateDerivatives must already be set by the caller (CtrlVar is a compile-time constant).
%
% Replaces the loop over integration points in uvhMatrixAssembly.m, ie the calls to uvhAssemblyIntPointImplicitSUPG_v2.m and
% uvhNodalLoopSSTREAM_v2.m, by a loop over elements, with the loop over integration points inside. For each element the nine
% nod x nod blocks of the element matrix and the element residuals are accumulated locally and then written to the outputs,
% which have the same layout as the arrays summed over integration points in uvhMatrixAssembly.m:
%
%   Tx, Fx, ... , Fh  : Nele x nod
%   Kxu, ... , Khh    : Nele x nod x nod, with K(e,j,i) the derivative of the residual of node i with respect to the unknown at node j
%   tauxInt, ...      : Nele x nip
%
% The terms are those of uvhAssemblyIntPointImplicitSUPG_v2.m and uvhNodalLoopSSTREAM_v2.m, written for one element and one
% integration point. The pointwise functions (EffectiveViscositySSTREAM, BasalDrag, SUPGtau, ...) are called with scalar arguments.
% Interpreted this is slow; it is compiled with MATLAB Coder (BuildUaCompiled.m).
%
%%

Nele=MUA.Nele; nod=MUA.nod; nip=MUA.nip; ndim=2;
assert(nod<=10) ; assert(nip<=37) ;   % upper bounds, so that the per-element arrays can be allocated statically


theta=CtrlVar.theta;
kH=CtrlVar.kH; Hh0=CtrlVar.Hh0; ThickMin=CtrlVar.ThickMin;
beta0=CtrlVar.SUPG.beta0; SpeedZero=CtrlVar.SpeedZero;
tauOption=CtrlVar.uvh.SUPG.tauCode;   % numeric code, see MakeKernelParamsFlex
tauMultiplier=CtrlVar.uvh.SUPG.tauMultiplier;

% options, constant over the loops
Melange=CtrlVar.IncludeMelangeModelPhysics;
UseLSF=CtrlVar.LevelSetMethod && CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback;
UsePenalty=isfield(CtrlVar,'ThicknessPenalty') && CtrlVar.ThicknessPenalty && ~(isfield(CtrlVar,'ThicknessPenaltyNodal') && CtrlVar.ThicknessPenaltyNodal);
UseBarrier=isfield(CtrlVar,'ThicknessBarrier') && isfield(CtrlVar,'ThicknessBarrierMassBalanceFeedbackCoeffLog') && CtrlVar.ThicknessBarrier;
hBarrier=0; pBarrier=0;   % defined on all paths (code generation)
if UseBarrier
    hBarrier=CtrlVar.ThickMin+eps(CtrlVar.ThickMin);
    pBarrier=abs(CtrlVar.ThicknessBarrierMassBalanceFeedbackCoeffLog);
end
OnlyDragAndViscosity=isfield(CtrlVar,'OnlyCalcBasalDragAndEffectiveViscosity') && CtrlVar.OnlyCalcBasalDragAndEffectiveViscosity;
ResetThickness=CtrlVar.ResetThicknessToMinThickness;
EnforcebAboveB=CtrlVar.Enforce_bAboveB;
hasq=~isempty(qv); hasmuk=~isempty(mukv); hasV0=~isempty(V0v);

% shape functions at the integration points, the same for all elements
N=zeros(nod,nip);
for Iint=1:nip
    N(:,Iint)=shape_fun(Iint,ndim,nod,MUA.points);
end
weights=MUA.weights;

Tx=zeros(Nele,nod); Fx=zeros(Nele,nod); Ty=zeros(Nele,nod); Fy=zeros(Nele,nod); Th=zeros(Nele,nod); Fh=zeros(Nele,nod);
% All element-matrix values in one array, element-major: Xmat(L,e), L = j + nod*(i-1) + nod^2*(c-1), with the nine blocks c in
% row-major order (u,v,h)x(u,v,h). The values are summed into the nonzeros of the sparse matrix with GatherSumCSC (see uvhCSCPattern.m).
if ~Ronly
    Xmat=zeros(9*nod*nod,Nele);
else
    Xmat=zeros(9*nod*nod,0);
end
tauxInt=zeros(Nele,nip); tauyInt=zeros(Nele,nip); etaInt=zeros(Nele,nip); HeInt=zeros(Nele,nip);

h0barr=0; h1barr=0; lambda_h=1;

parfor e=1:Nele
    idx=conn(e,:);   % nodes of element e (the nodal values are gathered here, step 3)
    % variables assigned only in some branches are initialised in each iteration (parfor)
    uoint=[]; voint=[]; Coint=[]; moint=[]; uaint=[]; vaint=[]; Caint=[]; maint=[];
    qint=[]; mukint=[]; V0int=[];
    bpn=zeros(1,nod);

    % nodal values of this element, 1 x nod
    hn=hv(idx).'; un=uv(idx).'; vn=vv(idx).'; bn=bv(idx).'; Bn=Bv(idx).'; Sn=Sv(idx).'; rhon=rhov(idx).';
    h0n=h0v(idx).'; u0n=u0v(idx).'; v0n=v0v(idx).';

    ConstrainedElement=any(hBCsMaskv(idx).'>0);   % element contains a node with h BCs or an active thickness constraint
    le=sqrt(2*MUA.EleAreas(e));

    % nodal derivative of b with respect to h (only needed for the Jacobian), see uvhAssemblyIntPointImplicitSUPG_v2.m
    if ~Ronly
        hrn=hn;
        if ResetThickness ; hrn(hn<ThickMin)=ThickMin ; end
        hfn=rhow*(Sn-Bn)./rhon;
        GFn=HeavisideApprox(kH,hrn-hfn,Hh0);
        deln=DiracDelta(kH,hrn-hfn,Hh0);
        bfn=Sn-rhon.*hrn/rhow;
        bpn=deln.*(Bn-bfn)-(1-GFn).*rhon/rhow;
        if EnforcebAboveB ; bpn(GFn.*Bn+(1-GFn).*bfn<Bn)=0 ; end
        if ResetThickness ; bpn(hn<ThickMin)=0 ; end
    end

    % element residuals and element matrices, accumulated over integration points. k(j,i) corresponds to K(e,j,i)
    tx=zeros(1,nod); fx=zeros(1,nod); ty=zeros(1,nod); fy=zeros(1,nod); th=zeros(1,nod); fh=zeros(1,nod);
    kxu=zeros(nod,nod); kxv=zeros(nod,nod); kxh=zeros(nod,nod);
    kyu=zeros(nod,nod); kyv=zeros(nod,nod); kyh=zeros(nod,nod);
    khu=zeros(nod,nod); khv=zeros(nod,nod); khh=zeros(nod,nod);

    for Iint=1:nip

        fun=N(:,Iint); funR=fun.';
        Dx=reshape(MUA.Deriv(e,1,:,Iint),1,nod);
        Dy=reshape(MUA.Deriv(e,2,:,Iint),1,nod);
        detJw=MUA.DetJ(e,Iint)*weights(Iint);

        % values at the integration point
        hint=hn*fun; uint=un*fun; vint=vn*fun;
        if Melange
            uoint=uov(idx).'*fun; voint=vov(idx).'*fun; uaint=uav(idx).'*fun; vaint=vav(idx).'*fun;
            Coint=Cov(idx).'*fun; moint=mov(idx).'*fun; Caint=Cav(idx).'*fun; maint=mav(idx).'*fun;
        end
        Cint=SmoothFloor(Cv(idx).'*fun,CtrlVar.Cmin,CtrlVar.CminWidth);
        mint=mv(idx).'*fun;
        if hasq ; qint=qv(idx).'*fun ; end
        if hasmuk ; mukint=mukv(idx).'*fun ; end
        if hasV0 ; V0int=V0v(idx).'*fun ; end
        AGlenint=SmoothFloor(AGlenv(idx).'*fun,CtrlVar.AGlenmin,CtrlVar.AGlenminWidth);
        nint=nv(idx).'*fun;

        h0int=h0n*fun; u0int=u0n*fun; v0int=v0n*fun;
        a1int=as1v(idx).'*fun+ab1v(idx).'*fun;
        a0int=as0v(idx).'*fun+ab0v(idx).'*fun;
        dadhint=0;

        if UseLSF
            LM=LSFMaskv(idx).'*fun;
            [abLSF,dadhLSF]=LevelSetMethodMassBalanceFeedback(CtrlVar,LM,hint);
            a1int=a1int+abLSF;
            dadhint=dadhint+dadhLSF;
        else
            LM=0;
        end

        if UsePenalty
            [aPenalty1,daPenaltydh1]=ThicknessPenaltyMassBalanceFeedbackFlex(CtrlVar,hint);
            aPenalty1=aPenalty1(1)*(~ConstrainedElement);      % (1): scalar, as hint is a scalar
            daPenaltydh1=daPenaltydh1(1)*(~ConstrainedElement);
            a1int=a1int+aPenalty1;
            dadhint=dadhint+daPenaltydh1;
        end

        if UseBarrier
            if hint > hBarrier && ~ConstrainedElement && LM < 0.1
                a1int=a1int-pBarrier*log(hint-hBarrier);
                dadhint=dadhint-pBarrier/(hint-hBarrier);
            end
        end

        Bint=Bn*fun; Sint=Sn*fun; rhoint=rhon*fun;
        Hint=Sint-Bint;
        hfint=rhow*Hint/rhoint;

        Heint=HeavisideApprox(kH,hint-hfint,Hh0);
        HEint=HeavisideApprox(kH,hfint-hint,Hh0);
        deltaint=DiracDelta(kH,hint-hfint,Hh0);
        Deltaint=DiracDelta(kH,hfint-hint,Hh0);
        Hposint=HeavisideApprox(kH,Hint,Hh0)*Hint;
        dint=HEint*rhoint*hint/rhow+Heint*Hposint;
        Dddhint=HEint*rhoint/rhow-Deltaint*hint*rhoint/rhow+deltaint*Hposint;

        % derivatives at the integration point
        dhdx=Dx*hn.'; dhdy=Dy*hn.';
        dbdx=Dx*bn.'; dbdy=Dy*bn.';
        dh0dx=Dx*h0n.'; dh0dy=Dy*h0n.';
        drhodx=Dx*rhon.'; drhody=Dy*rhon.';
        exx0=Dx*u0n.'; eyy0=Dy*v0n.';
        exx=Dx*un.'; eyy=Dy*vn.'; exy=0.5*(Dx*vn.'+Dy*un.');

        [etaint,Eint]=EffectiveViscositySSTREAM(CtrlVar,AGlenint,nint,exx,eyy,exy);
        [taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh]=...
            BasalDragFlex(CtrlVar,MUA,Heint,deltaint,hint,Bint,Hint,rhoint,rhow,uint,vint,Cint,mint,uoint,voint,Coint,moint,uaint,vaint,Caint,maint,qint,g,mukint,V0int);

        taux=taux(1); tauy=tauy(1); etaint=etaint(1);   % (1): scalars, as the inputs are scalars (sizes depend on run-time options in BasalDrag)
        tauxInt(e,Iint)=taux; tauyInt(e,Iint)=tauy; etaInt(e,Iint)=etaint; HeInt(e,Iint)=Heint;

        if OnlyDragAndViscosity
            continue
        end

        speed0=sqrt(u0int*u0int+v0int*v0int+SpeedZero^2);
        tau=SUPGtauFlex(CtrlVar,speed0,le,dt,tauOption,tauMultiplier);
        tau0=beta0*tau;

        % --- as uvhNodalLoopSSTREAM_v2.m, for one element ---
        qx1dx=rhoint*exx*hint+rhoint*uint*dhdx+drhodx*uint*hint;
        qy1dy=rhoint*eyy*hint+rhoint*vint*dhdy+drhody*vint*hint;
        qx0dx=rhoint*exx0*h0int+rhoint*u0int*dh0dx+drhodx*u0int*h0int;
        qy0dy=rhoint*eyy0*h0int+rhoint*v0int*dh0dy+drhody*v0int*h0int;

        SUPG=funR+0.5*tau0*(u0int*Dx+v0int*Dy);    % 1 x nod (i)

        if ~Ronly

            Eint=Eint(1); dtauxdu=dtauxdu(1); dtauxdv=dtauxdv(1); dtauydu=dtauydu(1); dtauydv=dtauydv(1); dtauxdh=dtauxdh(1); dtauydh=dtauydh(1);
            heW=hint*etaint*detJw;
            a1=hint*(4*exx+2*eyy);
            a2=2*hint*exy;
            a3=hint*(4*eyy+2*exx);

            Deu=Eint*((2*exx+eyy)*Dx+exy*Dy);      % 1 x nod (j)
            Dev=Eint*((2*eyy+exx)*Dy+exy*Dx);

            GuW=(a1*Dx+a2*Dy)*detJw;               % 1 x nod (i)
            GvW=(a3*Dy+a2*Dx)*detJw;

            tuuW=dtauxdu*detJw; tvvW=dtauydv*detJw; tuvW=dtauxdv*detJw; tvuW=dtauydu*detJw;

            t2coeff=ca*g*(rhoint*hint-rhow*dint*Dddhint);
            fxh=dtauxdh+ca*g*(rhoint-rhow*Dddhint)*dbdx-sa*g*rhoint;
            fyh=dtauydh+ca*g*(rhoint-rhow*Dddhint)*dbdy;

            PxW=((etaint*(4*exx+2*eyy)-t2coeff)*Dx+etaint*2*exy*Dy+fxh*funR)*detJw;   % 1 x nod (i)
            PyW=((etaint*(4*eyy+2*exx)-t2coeff)*Dy+etaint*2*exy*Dx+fyh*funR)*detJw;

            cDb=ca*g*(rhoint*hint-rhow*dint)*detJw;

            QuW=((rhoint*dhdx+drhodx*hint)*funR+rhoint*hint*Dx)*theta*detJw*dt;        % 1 x nod (j)
            QvW=((rhoint*dhdy+drhody*hint)*funR+rhoint*hint*Dy)*theta*detJw*dt;
            QhW=((rhoint-dt*theta*rhoint*dadhint-dt*theta*rhoint*dadhv(idx).'+dt*theta*rhoint*h1barr/lambda_h ...
                +dt*theta*(rhoint*exx+drhodx*uint+rhoint*eyy+drhody*vint)).*funR ...
                +dt*theta*rhoint*uint*Dx+dt*theta*rhoint*vint*Dy)*detJw;

            for i=1:nod
                dxi=Dx(i); dyi=Dy(i); PP=fun(i)*funR;
                kxu(:,i)=kxu(:,i)+(4*heW*dxi*Dx+heW*dyi*Dy+GuW(i)*Deu+tuuW*PP).';
                kyv(:,i)=kyv(:,i)+(4*heW*dyi*Dy+heW*dxi*Dx+GvW(i)*Dev+tvvW*PP).';
                kxv(:,i)=kxv(:,i)+(heW*(2*dxi*Dy+dyi*Dx)+GuW(i)*Dev+tuvW*PP).';
                kyu(:,i)=kyu(:,i)+(heW*(2*dyi*Dx+dxi*Dy)+GvW(i)*Deu+tvuW*PP).';
                kxh(:,i)=kxh(:,i)+(PxW(i)*funR).';
                kyh(:,i)=kyh(:,i)+(PyW(i)*funR).';
                kxh(:,i)=kxh(:,i)+((cDb*fun(i))*(Dx.*bpn)).';
                kyh(:,i)=kyh(:,i)+((cDb*fun(i))*(Dy.*bpn)).';
                khu(:,i)=khu(:,i)+(SUPG(i)*QuW).';
                khv(:,i)=khv(:,i)+(SUPG(i)*QvW).';
                khh(:,i)=khh(:,i)+(SUPG(i)*QhW).';
            end

        end

        % residuals
        c1x=(-ca*g*(rhoint*hint-rhow*dint)*dbdx+rhoint*g*hint*sa)*detJw;
        c1y=(-ca*g*(rhoint*hint-rhow*dint)*dbdy)*detJw;
        c2=(0.5*ca*g*(rhoint*hint^2-rhow*dint^2))*detJw;
        c3x=(hint*etaint*(4*exx+2*eyy))*detJw;
        c3y=(hint*etaint*(4*eyy+2*exx))*detJw;
        c4=(hint*etaint*2*exy)*detJw;
        txW=taux*detJw;
        tyW=tauy*detJw;

        tx=tx+c3x*Dx+c4*Dy+txW*funR;
        fx=fx+c1x*funR+c2*Dx;
        ty=ty+c3y*Dy+c4*Dx+tyW*funR;
        fy=fy+c1y*funR+c2*Dy;

        qterm=dt*(theta*qx1dx+(1-theta)*qx0dx+theta*qy1dy+(1-theta)*qy0dy);
        dhdtq=rhoint*(h0int-hint+dt*(1-theta)*h0barr+dt*theta*h1barr);
        accterm=dt*rhoint*((1-theta)*a0int+theta*a1int);

        th=th+(-dhdtq*detJw)*SUPG;
        fh=fh+((accterm-qterm)*detJw)*SUPG;

    end

    Tx(e,:)=tx; Fx(e,:)=fx; Ty(e,:)=ty; Fy(e,:)=fy; Th(e,:)=th; Fh(e,:)=fh;
    if ~Ronly
        Xmat(:,e)=[kxu(:);kxv(:);kxh(:);kyu(:);kyv(:);kyh(:);khu(:);khv(:);khh(:)];   % element-major (contiguous for each element)
    end

end

end
