function [Tx,Fx,Ty,Fy,vs]=uvElementKernel(MUA,Deriv,DetJ,conn,...
    sv,hv,uv,vv,AGlenv,nv,Cv,mv,qv,mukv,V0v,Bv,Sv,rhov,...
    uov,vov,Cov,mov,uav,vav,Cav,mav,...
    CtrlVar,rhow,g,Ronly,ca,sa,invPtr,invIdx,perm)

%% Element-based uv assembly kernel, compiled with MATLAB Coder (30 Oct 2026); see README.md in this folder
%
% Called by uvMatrixAssemblyCompiled.m through the MEX files built by BuildUaCompiled.m (entry points uvKernelFull.m, uvKernelRonly.m).
% Replaces the loop over integration points in uvMatrixAssemblySSTREAM.m by a loop over elements (parfor, OpenMP), with the loop over
% integration points inside. The terms are those of uvMatrixAssemblySSTREAM.m, written for one element and one integration point; the
% pointwise functions (BasalDrag, EffectiveViscositySSTREAM, SmoothFloor, HeavisideApprox, DiracDelta) are called with scalar arguments.
%
% Inputs:
%   MUA        : Nele, nod, nip, points, weights
%   Deriv, DetJ: MUA.Deriv and MUA.DetJ, passed as separate arrays (fields of a structure input are copied into a MEX file, arrays are not)
%   conn       : MUA.connectivity; the nodal values of each element are gathered from the nodal vectors (Nnodes x 1) here
%   CtrlVar    : run-time options and parameters (MakeKernelParamsFlex.m); CtrlVar.EffectiveViscosity.CalculateDerivatives is set by the
%                entry points
%   Ronly      : compile-time constant (residual only)
%   invPtr, invIdx : from CompiledAssemblyPattern.m, for summing the element-matrix values into the nonzeros
%   perm       : from the assembly pattern (MUA.uvAssemblyPattern.perm), for the symmetrisation of the nonzeros
%
% Outputs:
%   Tx, Fx, Ty, Fy : element residuals (internal and external forces), Nele x nod
%   vs             : the nonzero values of the tangent matrix, in CSC order, symmetrised as in uvMatrixAssemblySSTREAM.m (empty if Ronly).
%                    The element-matrix values are kept internally (Xmat, element-major) and summed into the nonzeros in parallel over the
%                    nonzeros, in the order given by invIdx
%
% Outputs that are completely overwritten are not zero-filled (coder.nullcopy). Interpreted, this function is slow.
%
%%

Nele=MUA.Nele; nod=MUA.nod; nip=MUA.nip; ndim=2;
assert(nod<=10) ; assert(nip<=37) ;   % upper bounds, so that the per-element arrays can be allocated statically

kH=CtrlVar.kH; Hh0=CtrlVar.Hh0;

% options, constant over the loops
Melange=CtrlVar.IncludeMelangeModelPhysics;
DragDerivatives=double(~Ronly);   % BasalDrag: derivatives with respect to u and v, unless Ronly
hasq=~isempty(qv); hasmuk=~isempty(mukv); hasV0=~isempty(V0v);

% shape functions at the integration points, the same for all elements
N=zeros(nod,nip);
for Iint=1:nip
    N(:,Iint)=shape_fun(Iint,ndim,nod,MUA.points);
end
weights=MUA.weights;

% (outputs that are completely overwritten are not zero-filled)
Tx=coder.nullcopy(zeros(Nele,nod)); Fx=coder.nullcopy(zeros(Nele,nod)); Ty=coder.nullcopy(zeros(Nele,nod)); Fy=coder.nullcopy(zeros(Nele,nod));
% All element-matrix values in one array, element-major: Xmat(L,e), L = j + nod*(i-1) + nod^2*(c-1), with the four blocks c in the
% order of uvMatrixAssemblySSTREAM.m: uu, vv, uv, vu (Xval=[d1d1(:);d2d2(:);d1d2(:);d2d1(:)]), on which the assembly pattern is based.
if ~Ronly
    Xmat=coder.nullcopy(zeros(4*nod*nod,Nele));
else
    Xmat=zeros(4*nod*nod,0);
end

parfor e=1:Nele
    idx=conn(e,:);   % nodes of element e
    % variables assigned only in some branches are initialised in each iteration (parfor)
    uoint=0; voint=0; Coint=0; moint=0; uaint=0; vaint=0; Caint=0; maint=0;   % used only with melange
    qint=NaN; mukint=NaN; V0int=NaN;   % used only by the sliding laws that need them

    % nodal values of this element, 1 x nod
    hn=hv(idx).'; sn=sv(idx).'; un=uv(idx).'; vn=vv(idx).'; Bn=Bv(idx).'; Sn=Sv(idx).'; rhon=rhov(idx).';

    % element residuals and element matrices, accumulated over integration points. k(j,i) corresponds to d(e,j,i) in
    % uvMatrixAssemblySSTREAM.m (i: test function, j: unknown)
    tx=zeros(1,nod); fx=zeros(1,nod); ty=zeros(1,nod); fy=zeros(1,nod);
    kuu=zeros(nod,nod); kvv=zeros(nod,nod); kuv=zeros(nod,nod); kvu=zeros(nod,nod);

    for Iint=1:nip

        fun=N(:,Iint); funR=fun.';
        Dx=reshape(Deriv(e,1,:,Iint),1,nod);
        Dy=reshape(Deriv(e,2,:,Iint),1,nod);
        detJw=DetJ(e,Iint)*weights(Iint);

        % values at the integration point
        uint=un*fun; vint=vn*fun;
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

        Bint=Bn*fun; Sint=Sn*fun;
        Hint=Sint-Bint;
        hint=hn*fun;
        rhoint=rhon*fun;

        % flotation thickness, Heaviside and Dirac functions, and d, at the integration point (as in uvMatrixAssemblySSTREAM.m)
        hfint=rhow*Hint/rhoint;
        Heint=HeavisideApprox(kH,hint-hfint,Hh0);
        HEint=HeavisideApprox(kH,hfint-hint,Hh0);
        deltaint=DiracDelta(kH,hint-hfint,Hh0);
        Hposint=HeavisideApprox(kH,Hint,Hh0)*Hint;
        dint=HEint*rhoint*hint/rhow+Heint*Hposint;

        % derivatives at the integration point
        dsdx=Dx*sn.'; dsdy=Dy*sn.';
        dhdx=Dx*hn.'; dhdy=Dy*hn.';
        exx=Dx*un.'; eyy=Dy*vn.'; exy=0.5*(Dx*vn.'+Dy*un.');

        [taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv]=...
            BasalDrag(CtrlVar,DragDerivatives,Heint,deltaint,hint,Hint,rhoint,rhow,g,uint,vint,Cint,mint,qint,mukint,V0int,uoint,voint,Coint,moint,uaint,vaint,Caint,maint);
        [etaint,Eint]=EffectiveViscositySSTREAM(CtrlVar,AGlenint,nint,exx,eyy,exy);
        taux=taux(1); tauy=tauy(1); etaint=etaint(1);   % (1): scalars, as the inputs are scalars

        dbdx=dsdx-dhdx;
        dbdy=dsdy-dhdy;

        % tangent matrix
        if ~Ronly
            Eint=Eint(1); dtauxdu=dtauxdu(1); dtauxdv=dtauxdv(1); dtauydu=dtauydu(1); dtauydv=dtauydv(1);
            he=hint*etaint*detJw;
            tuu=dtauxdu*detJw; tvv=dtauydv*detJw; tuv=dtauxdv*detJw; tvu=dtauydu*detJw;
            a1=hint*(4*exx+2*eyy);
            a2=2*hint*exy;
            a3=hint*(4*eyy+2*exx);
            Deu=Eint*((2*exx+eyy)*Dx+exy*Dy);      % 1 x nod (j)
            Dev=Eint*((2*eyy+exx)*Dy+exy*Dx);
            GuW=(a1*Dx+a2*Dy)*detJw;               % 1 x nod (i)
            GvW=(a3*Dy+a2*Dx)*detJw;
            for i=1:nod
                dxi=Dx(i); dyi=Dy(i); PP=fun(i)*funR;
                kuu(:,i)=kuu(:,i)+(4*he*dxi*Dx+he*dyi*Dy+tuu*PP+GuW(i)*Deu).';
                kvv(:,i)=kvv(:,i)+(4*he*dyi*Dy+he*dxi*Dx+tvv*PP+GvW(i)*Dev).';
                kuv(:,i)=kuv(:,i)+(he*(2*dxi*Dy+dyi*Dx)+tuv*PP+GuW(i)*Dev).';
                kvu(:,i)=kvu(:,i)+(he*(2*dyi*Dx+dxi*Dy)+tvu*PP+GvW(i)*Deu).';
            end
        end

        % residual (internal forces T) and external forces (F)
        qx =(-g*(rhoint*hint-rhow*dint)*dbdx*ca+rhoint*g*hint*sa)*detJw;
        qy =(-g*(rhoint*hint-rhow*dint)*dbdy*ca)*detJw;
        p2 =(0.5*g*ca*(rhoint*hint^2-rhow*dint^2))*detJw;
        p3x=(hint*etaint*(4*exx+2*eyy))*detJw;
        p3y=(hint*etaint*(4*eyy+2*exx))*detJw;
        p4 =(hint*etaint*2*exy)*detJw;
        txW=taux*detJw;
        tyW=tauy*detJw;

        tx=tx+p3x*Dx+p4*Dy+txW*funR;
        fx=fx+qx*funR+p2*Dx;
        ty=ty+p3y*Dy+p4*Dx+tyW*funR;
        fy=fy+qy*funR+p2*Dy;

    end

    Tx(e,:)=tx; Fx(e,:)=fx; Ty(e,:)=ty; Fy(e,:)=fy;
    if ~Ronly
        Xmat(:,e)=[kuu(:);kvv(:);kuv(:);kvu(:)];   % element-major (contiguous for each element), blocks as Xval in uvMatrixAssemblySSTREAM.m
    end

end


% Sum the element-matrix values into the nonzeros of the sparse matrix, in parallel over the nonzeros (each written by one thread),
% in the order given by invIdx (see CompiledAssemblyPattern.m), and symmetrise: vs=0.5*(vs+vs(perm)), as in uvMatrixAssemblySSTREAM.m.
% The element-matrix values Xmat stay inside the MEX file.
if ~Ronly
    nk=numel(invPtr)-1;
    v0=coder.nullcopy(zeros(nk,1));
    parfor k=1:nk
        s=0;
        for p=invPtr(k)+1:invPtr(k+1)
            s=s+Xmat(invIdx(p));
        end
        v0(k)=s;
    end
    vs=coder.nullcopy(zeros(nk,1));
    parfor k=1:nk
        vs(k)=0.5*(v0(k)+v0(perm(k)));
    end
else
    vs=zeros(0,1);
end
end
