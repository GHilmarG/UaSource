




function [UserVar,RunInfo,R,KdFuvhduvh,tauxInt,tauyInt,etaInt,HeInt]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)

% [UserVar,RunInfo,R,K,Tint,Fext]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,ZeroFields)
%
%
% Does not depend on s or b, only h, S and B
%
% K is dFuvh/duvh
%
%

%
% Tint   : vector of internal nodal forces
% Fext   : vector of external nodal forces


narginchk(8,8)
nargoutchk(4,8)



%
% K=
%    [Kxu Kxv Kxh]
%    [Kyu Kyv Kyh]
%    [Khu Khv Khh]

ZeroFields=CtrlVar.uvhMatrixAssembly.ZeroFields;
Ronly=CtrlVar.uvhMatrixAssembly.Ronly;

if Ronly
    KdFuvhduvh=[];
end


if Ronly
    CtrlVar.BasalDrag.CalculateDerivatives=false;  
    CtrlVar.EffectiveViscosity.CalculateDerivatives=false;
else
    CtrlVar.BasalDrag.CalculateDerivatives=true;
    CtrlVar.EffectiveViscosity.CalculateDerivatives=true;

end
%
% if nargin<7
%     ZeroFields=false;
% end
%
%
% if nargout==3
%     Ronly=1;
% else
%     Ronly=0;
% end


if ZeroFields

    % I'm using this to come up with a reasonable normalizing factor to the residuals.  The uv side of things is clear and there
    % I set u=v=0 and this ensures that all 'internal' nodal forces are zero. The h side is less clear. and I've struggled with
    % finding a sensible normalizing factor for this term. If I set u=v=0 to get the sensible uv normalization, then dqdx=0. I
    % can have a situation where a=0 and if h1=h0 then the initial estimate for dh/dt=0. So all terms are then zero.
    %
    % One approach is to create a scale for dh/dt by using ThickMin/dt , or alternatively just make this term numerically equal
    % to 1, i.e. no normalization apart in the norm. This can be achieved by setting the surface mass balance to 1.
    %
    % Possibly best to normalize with the prescribed mass balance. The justification for doing that would be that the mass
    % balance is externally prescribed to the mass conservation equation, similar to how the body force is external load in the
    % momentum equation. The h-equation is then sufficiently accurately solved once the ratio:
    %
    %    (dh/dt-dq/dx-a)/a
    %
    % is small. However the issue is that this will not work if a=0 everywhere.
    %
    % I therefore somewhat arbitrarily add 1 to this
    %
    %     (dh/dt-dq/dx-a)/(abs(a)+1)
    %
    % Note that I can use the abs here because this is only used for the normalization factor, which is later squired.
    %

    F1.ub=F1.ub*0; F1.vb=F1.vb*0;
    F0.ub=F0.ub*0; F0.vb=F0.vb*0;
    %
    % How to normalize the mass conservation term?
    %
    % Idea1) set a=1 as a normalizing factor
    % The issue with this is that the accterm -> 0 as dt->0
    % because
    % accterm=  dt*rhoint.*((1-theta)*a0int+theta*a1int).*SUPG;
    % so the normalisation factor goes to zero with dt


    F1.h=F0.h;  % this leads to a dh/dt=0 at the beginning
    F1.as=abs(F1.ab)+1;  F1.ab=abs(F1.ab);  % I can use abs here because this is just for the normalization factor which is squared.
    F0.as=abs(F0.ab)+1;  F0.ab=abs(F0.ab);



    % I can solve this by dividing with dt again as I calculate the normalization factor in the const function. This means that
    % the normalization is independent of dt
    %
    % Possibly it would be better to solve directly for dh/dt, then at least the units of the rhs are identical for all unknowns
    %
    % On the other hand this can hardly be too much of an issue as the dh/dt equation is linear in h and all the residuals will
    % be caused by the u v residuals.
    %


end



if ~CtrlVar.IncludeMelangeModelPhysics
    uonod=[];
    vonod=[];
    Conod=[];
    monod=[];
    uanod=[];
    vanod=[];
    Canod=[];
    manod=[];
end




if any(isnan(F1.ub))
    fprintf(CtrlVar.fidlog,' NaN in u on input to uvhMatrixAssembly \n');
end
if any(isnan(F1.vb)) ;  fprintf(CtrlVar.fidlog,' NaN in v on input to uvhMatrixAssembly \n'); end
if any(isnan(F1.h)) ;  fprintf(CtrlVar.fidlog,' NaN in h on input to uvhMatrixAssembly \n'); end
if any(isnan(F0.ub)) ;  fprintf(CtrlVar.fidlog,' NaN in u0 on input to uvhMatrixAssembly \n'); end
if any(isnan(F0.vb)) ;  fprintf(CtrlVar.fidlog,' NaN in v0 on input to uvhMatrixAssembly \n'); end
if any(isnan(F0.h)) ;  fprintf(CtrlVar.fidlog,' NaN in h0 on input to uvhMatrixAssembly \n'); end

g=F1.g ;
alpha=F1.alpha;
rhow=F1.rhow;
dt=CtrlVar.dt;

temp=CtrlVar.ResetThicknessToMinThickness;
if ~CtrlVar.ResetThicknessInNonLinLoop
    CtrlVar.ResetThicknessToMinThickness=0;
end


[F1.b,F1.s]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);   % don't update h outside of loop, just get new values for s and b

CtrlVar.ResetThicknessToMinThickness=temp;

if CtrlVar.MassBalanceGeometryFeedback>=2  && ~ZeroFields

    rdamp=CtrlVar.MassBalanceGeometryFeedbackDamping;
    if rdamp~=0
        as1Old=F1.as ; ab1Old=F1.ab;
    end
    CtrlVar.time=CtrlVar.time+CtrlVar.dt;
    [UserVar,F1]=GetMassBalance(UserVar,CtrlVar,MUA,F1);
    CtrlVar.time=CtrlVar.time-CtrlVar.dt;
    switch CtrlVar.MassBalanceGeometryFeedback

        case 2
            dadh=zeros(MUA.Nnodes,1);
        case 3
            dadh=F1.dasdh+F1.dabdh;
    end


    if rdamp~=0
        % I don't account for a potential dependency of as and ab
        % on h in the Hessian, so may need to dampen these changes
        F1.as=(1-rdamp)*F1.as+rdamp*as1Old;
        F1.ab=(1-rdamp)*F1.ab+rdamp*ab1Old;
    end
else
    dadh=zeros(MUA.Nnodes,1);
end



if CtrlVar.LevelSetMethod  &&  CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback  && ~isempty(F1.LSF)
    if isempty(F1.LSFMask)
        F1.LSFMask=CalcMeshMask(CtrlVar,MUA,F1.LSF,0);

    end
    LSFMask=F1.LSFMask.NodesOut ; % This is the 'strictly' definition
else

    LSFMask=zeros(MUA.Nnodes,1) ;
end


ndim=2;  neq=3*MUA.Nnodes;
neqx=MUA.Nnodes ;

hnod=reshape(F1.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(F1.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F1.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

LSFMasknod=reshape(LSFMask(MUA.connectivity,1),MUA.Nele,MUA.nod);

hBCsMasknod=zeros(MUA.Nnodes,1) ;
hBCsMasknod(BCs1.hFixedNode)=1;
hBCsMasknod(BCs1.hPosNode)=1;
hBCsMasknod=reshape(hBCsMasknod(MUA.connectivity,1),MUA.Nele,MUA.nod);

% hBCs=BCs.hFixedNode;BCs.hPosNode]

if CtrlVar.IncludeMelangeModelPhysics

    uonod=reshape(F1.uo(MUA.connectivity,1),MUA.Nele,MUA.nod);
    vonod=reshape(F1.vo(MUA.connectivity,1),MUA.Nele,MUA.nod);

    uanod=reshape(F1.ua(MUA.connectivity,1),MUA.Nele,MUA.nod);
    vanod=reshape(F1.va(MUA.connectivity,1),MUA.Nele,MUA.nod);

end



Cnod=reshape(F1.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F1.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

if ~isempty(F1.q)
    qnod=reshape(F1.q(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    qnod=[];
end

if ~isempty(F1.muk)
    muknod=reshape(F1.muk(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    muknod=[];
end

if ~isempty(F1.V0)
    V0nod=reshape(F1.V0(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    V0nod=[];
end


if CtrlVar.IncludeMelangeModelPhysics

    Conod=reshape(F1.Co(MUA.connectivity,1),MUA.Nele,MUA.nod);
    monod=reshape(F1.mo(MUA.connectivity,1),MUA.Nele,MUA.nod);
    Canod=reshape(F1.Ca(MUA.connectivity,1),MUA.Nele,MUA.nod);
    manod=reshape(F1.ma(MUA.connectivity,1),MUA.Nele,MUA.nod);

end

AGlennod=reshape(F1.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F1.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

Snod=reshape(F1.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(F1.B(MUA.connectivity,1),MUA.Nele,MUA.nod);

h0nod=reshape(F0.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
u0nod=reshape(F0.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
v0nod=reshape(F0.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
as0nod=reshape(F0.as(MUA.connectivity,1),MUA.Nele,MUA.nod);
ab0nod=reshape(F0.ab(MUA.connectivity,1),MUA.Nele,MUA.nod);

as1nod=reshape(F1.as(MUA.connectivity,1),MUA.Nele,MUA.nod);
ab1nod=reshape(F1.ab(MUA.connectivity,1),MUA.Nele,MUA.nod);
dadhnod=reshape(dadh(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F1.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
bnod=reshape(F1.b(MUA.connectivity,1),MUA.Nele,MUA.nod);
%dudtnod=reshape(F1.dubdt(MUA.connectivity,1),MUA.Nele,MUA.nod);
%dvdtnod=reshape(F1.dvbdt(MUA.connectivity,1),MUA.Nele,MUA.nod);


ca=cos(alpha); sa=sin(alpha);



if ~Ronly
    Kxu=zeros(MUA.Nele,MUA.nod,MUA.nod); Kxv=zeros(MUA.Nele,MUA.nod,MUA.nod);  Kxh=zeros(MUA.Nele,MUA.nod,MUA.nod);
    Kyu=zeros(MUA.Nele,MUA.nod,MUA.nod); Kyv=zeros(MUA.Nele,MUA.nod,MUA.nod);  Kyh=zeros(MUA.Nele,MUA.nod,MUA.nod);
    Khu=zeros(MUA.Nele,MUA.nod,MUA.nod); Khv=zeros(MUA.Nele,MUA.nod,MUA.nod);  Khh=zeros(MUA.Nele,MUA.nod,MUA.nod);
else
    Kxu=[]; Kxv=[];  Kxh=[];
    Kyu=[]; Kyv=[];  Kyh=[];
    Khu=[]; Khv=[];  Khh=[];
end

Tx=zeros(MUA.Nele,MUA.nod);  Ty=zeros(MUA.Nele,MUA.nod); Fx=zeros(MUA.Nele,MUA.nod);  Fy=zeros(MUA.Nele,MUA.nod); Th=zeros(MUA.Nele,MUA.nod);  Fh=zeros(MUA.Nele,MUA.nod);

Tx0=Tx ;Fx0=Fx; Ty0=Ty ; Fy0=Fy; Th0=Th ;Fh0=Fh;
Kxu0=Kxu ; Kxv0=Kxv ; Kyu0=Kyu ; Kyv0=Kyv ; Kxh0=Kxh ; Kyh0=Kyh ; Khu0=Khu ; Khv0=Khv ; Khh0=Khh;

if nargout> 4
    tauxInt=zeros(MUA.Nele,MUA.nip) ;
    tauyInt=zeros(MUA.Nele,MUA.nip) ;
    etaInt=zeros(MUA.Nele,MUA.nip) ;
    HeInt=zeros(MUA.Nele,MUA.nip) ;
else
    tauxInt=[];
    tauyInt=[];
    etaInt=[];
    HeInt=[];
end


% adding contribution from each form function for every element. The number of form functions equals the number of nodes, so
% the outputs have the dimension ele x nodes
for Iint=1:MUA.nip

    [Tx1,Fx1,Ty1,Fy1,Th1,Fh1,Kxu1,Kxv1,Kyu1,Kyv1,Kxh1,Kyh1,Khu1,Khv1,Khh1,tauxI,tauyI,etaI,HeI]=...
        uvhAssemblyIntPointImplicitSUPG(Iint,ndim,MUA,...
        bnod,hnod,unod,vnod,AGlennod,nnod,Cnod,mnod,qnod,muknod,V0nod,h0nod,u0nod,v0nod,as0nod,ab0nod,as1nod,ab1nod,dadhnod,Bnod,Snod,rhonod,...
        LSFMasknod,hBCsMasknod,...
        uonod,vonod,Conod,monod,uanod,vanod,Canod,manod,...
        CtrlVar,rhow,g,Ronly,ca,sa,dt,...
        Tx0,Fx0,Ty0,Fy0,Th0,Fh0,Kxu0,Kxv0,Kyu0,Kyv0,Kxh0,Kyh0,Khu0,Khv0,Khh0);

    Tx=Tx+Tx1;  Fx=Fx+Fx1;
    Ty=Ty+Ty1;  Fy=Fy+Fy1;
    Th=Th+Th1;  Fh=Fh+Fh1;

    Kxu=Kxu+Kxu1;        Kxv=Kxv+Kxv1;
    Kyu=Kyu+Kyu1;        Kyv=Kyv+Kyv1;
    Kxh=Kxh+Kxh1;        Kyh=Kyh+Kyh1;
    Khu=Khu+Khu1;        Khv=Khv+Khv1;        Khh=Khh+Khh1;

    if nargout> 4
        tauxInt(:,Iint)=tauxI;
        tauyInt(:,Iint)=tauyI;
        etaInt(:,Iint)=etaI;
        HeInt(:,Iint)=HeI;
    end

end


if ~isfield(CtrlVar,"OnlyCalcBasalDragAndEffectiveViscosity")
    CtrlVar.OnlyCalcBasalDragAndEffectiveViscosity=false;
end


if CtrlVar.OnlyCalcBasalDragAndEffectiveViscosity

    R=[] ; KdFuvhduvh=[] ;
    return

end




%% assemble right-hand side (fewer sparse calls, 30 Jan 2023)

iR=zeros(MUA.nod*MUA.Nele*3,1,"uint32");

Tval=zeros(MUA.nod*MUA.Nele*3,1);
Fval=zeros(MUA.nod*MUA.Nele*3,1);
istak=0;

for Inod=1:MUA.nod


    iR(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);


    Tval(istak+1:istak+MUA.Nele)=Tx(:,Inod);
    Fval(istak+1:istak+MUA.Nele)=Fx(:,Inod);

    istak=istak+MUA.Nele;
    iR(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx;
    Tval(istak+1:istak+MUA.Nele)=Ty(:,Inod);
    Fval(istak+1:istak+MUA.Nele)=Fy(:,Inod);

    istak=istak+MUA.Nele;
    iR(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx;
    Tval(istak+1:istak+MUA.Nele)=Th(:,Inod);
    Fval(istak+1:istak+MUA.Nele)=Fh(:,Inod);

    istak=istak+MUA.Nele;


end


Tint=accumarray(iR,Tval,[neq 1]);   % full: was sparse with density ~0.92
Fext=accumarray(iR,Fval,[neq 1]);   % full: was sparse with density ~0.92




%%

R=Tint-Fext;
if any(isnan(R)) ||  any(isnan(Tint)) || any(isnan(Fext))

    fprintf("nan in R or Tint or Fext")

end
% R=Tint-Fext;
% Tint=[Tx ; Ty ; Th] ;
% Rint=[Fx ; Fy ; Fh] ;

if ~Ronly




    if isfield(MUA,"uvhAssemblyPattern") && ~isempty(MUA.uvhAssemblyPattern)
        % cached sparsity pattern.  Block order must match
        % MUA.uvhAssemblyPattern.blocks, i.e. row-major (u,v,h)x(u,v,h).
        P=MUA.uvhAssemblyPattern;
        Xval=[Kxu(:);Kxv(:);Kxh(:);Kyu(:);Kyv(:);Kyh(:);Khu(:);Khv(:);Khh(:)];

        KdFuvhduvh=sparse(P.i0,P.j0,accumarray(P.map,Xval,[P.nk 1]),neq,neq);
    else
        %Iind=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1); Jind=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1);
        Iind=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1,"uint32"); Jind=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1,"uint32");

        Xval=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1);
        istak=0;
        for Inod=1:MUA.nod
            for Jnod=1:MUA.nod

                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Kxu(:,Jnod,Inod);
                istak=istak+MUA.Nele;

                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Kxv(:,Jnod,Inod);
                istak=istak+MUA.Nele;

                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Kxh(:,Jnod,Inod);
                istak=istak+MUA.Nele;

                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Kyu(:,Jnod,Inod);
                istak=istak+MUA.Nele;

                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Kyv(:,Jnod,Inod);
                istak=istak+MUA.Nele;

                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Kyh(:,Jnod,Inod);
                istak=istak+MUA.Nele;

                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Khu(:,Jnod,Inod);
                istak=istak+MUA.Nele;

                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Khv(:,Jnod,Inod);
                istak=istak+MUA.Nele;

                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Khh(:,Jnod,Inod);
                istak=istak+MUA.Nele;
            end

        end


        KdFuvhduvh=sparse(Iind,Jind,Xval,neq,neq);
    end



end



minh=min(F1.h);


if minh<2*CtrlVar.ThickMin && CtrlVar.InfoLevelNonLinIt>1000   % if min thickness is approaching ThickMin give some information on h within NR loop
    msg=sprintf('In NRuvh loop, assembly stage: min(h) %-f \t max(h) %-g \n ',minh,max(F1.h)) ;
    fprintf(CtrlVar.fidlog,msg) ;
end

if ~Ronly
    if full(any(isnan(diag(KdFuvhduvh))))
        error(' NaN in K ' ) ;
    end
end

if any(isnan(R))

    error("uvhMatrixAssembly:NaNinR"," NaN in R " ) ;
end
end







