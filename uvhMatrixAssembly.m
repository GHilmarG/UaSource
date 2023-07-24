function [UserVar,RunInfo,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1)

% [UserVar,RunInfo,R,K,Tint,Fext]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,ZeroFields)
%
%
% Does not depend on s or b, only h, S and B


%
% Tint   : vector of internal nodal forces
% Fext   : vector of external nodal forces


narginchk(6,6)
nargoutchk(4,4)

%
% K= [Kxu Kxv Kxh]
%    [Kyu Kyv Kyh]
%    [Khu Khv Khh]

ZeroFields=CtrlVar.uvhMatrixAssembly.ZeroFields;
Ronly=CtrlVar.uvhMatrixAssembly.Ronly;

if Ronly
    K=[];
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
    
    % I'm using this to come up with a reasonable normalizing factor to the
    % residuals.  The uv side of things is clear and there I set u=v=0 and this
    % ensures that all 'interanal' nodal forces are zero. The h side is less clear.
    % and I've stuggled with finding a sensible normalising factor for this term.
    % If I set u=v=0 to get the sensible uv normalisation, then dqdx=0. I can have a
    % situation where a=0 and if h1=h0 then the initial estimate for dh/dt=0. So all
    % terms are then zero. 
    %
    % One approach is to create a scale for dh/dt by using ThickMin/dt , or
    % alternativily just make this term numerically equal to 1, i.e. no
    % normalisatiion apart in the norm. This can be achived by setting the surface
    % mass balance to 1.
    
    F1.ub=F1.ub*0; F1.vb=F1.vb*0;
    F0.ub=F0.ub*0; F0.vb=F0.vb*0;  
    
    % How to normalize the mass conservation term?
    %
    % Idea1) set a=1 as a normalizing factor
    % The issue with this is that the accterm -> 0 as dt->0 
    % because
    % accterm=  dt*rhoint.*((1-theta)*a0int+theta*a1int).*SUPG;
    % so the normalisation factor goes to zero with dt
    F1.h=F0.h;  % this leads to a dh/dt=0 at the beginning
    F1.as=F1.ab*0+1;  F1.ab=F1.ab*0;  
    F0.as=F0.ab*0+1;  F0.ab=F0.ab*0;
    
    % I can solve this by dividing with dt again as I calculate the normalisation factor in the const
    % function. This means that the normalisation is independent of dt
    % 
    % Possibly it would be better to solve directly for dh/dt, then at
    % least the units of the rhs are identical for all unknonws
    %
    % On the other hand this can hardly be too much of an issue as the
    % dh/dt equation is linear in h and all the residuals will be caused by
    % the u v residuals.
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




if any(isnan(F1.ub)) ;  fprintf(CtrlVar.fidlog,' NaN in u on input to uvhMatrixAssembly \n'); end
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
    LSFMask=F1.LSFMask.NodesOut ; % This is the 'strickly' definition
else
    LSFMask=zeros(MUA.Nnodes,1) ;
end


ndim=2;  neq=3*MUA.Nnodes;
neqx=MUA.Nnodes ;

hnod=reshape(F1.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(F1.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F1.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

LSFMasknod=reshape(LSFMask(MUA.connectivity,1),MUA.Nele,MUA.nod);



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


if CtrlVar.uvGroupAssembly

    hfnod=rhow*(Snod-Bnod)./rhonod;

    deltanod=DiracDelta(CtrlVar.kH,hnod-hfnod,CtrlVar.Hh0);
    Deltanod=DiracDelta(CtrlVar.kH,hfnod-hnod,CtrlVar.Hh0);

    Henod = HeavisideApprox(CtrlVar.kH,hnod-hfnod,CtrlVar.Hh0);
    HEnod = HeavisideApprox(CtrlVar.kH,hfnod-hnod,CtrlVar.Hh0);
    Hnod=Snod-Bnod; 
    Hposnod = HeavisideApprox(CtrlVar.kH,Hnod,CtrlVar.Hh0).*Hnod;

    %    dnod = Hposnod.*(Snod-bnod);  % draft

    dnod=HEnod.*rhonod.*hnod/rhow+Henod.*Hposnod ;  % definition of d
    Dddhnod=HEnod.*rhonod/rhow-Deltanod.*hnod.*rhonod/rhow+deltanod.*Hposnod; % derivative of dnod with respect to hnod

else


    Henod=[] ;  deltanod=[] ;   Hposnod=[] ;  dnod=[];   Dddhnod=[];

end






% [points,weights]=sample('triangle',nip,ndim);

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


if CtrlVar.Parallel.uvhAssembly.parfor.isOn
    
    parfor Iint=1:MUA.nip
        
        
        [Tx1,Fx1,Ty1,Fy1,Th1,Fh1,Kxu1,Kxv1,Kyu1,Kyv1,Kxh1,Kyh1,Khu1,Khv1,Khh1]=...
            uvhAssemblyIntPointImplicitSUPG(Iint,ndim,MUA,...
            bnod,hnod,unod,vnod,AGlennod,nnod,Cnod,mnod,qnod,muknod,h0nod,u0nod,v0nod,as0nod,ab0nod,as1nod,ab1nod,dadhnod,Bnod,Snod,rhonod,...
            Henod,deltanod,Hposnod,dnod,Dddhnod,...
            LSFMasknod,...
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
        
    end
    
else
    
    for Iint=1:MUA.nip
        
        [Tx1,Fx1,Ty1,Fy1,Th1,Fh1,Kxu1,Kxv1,Kyu1,Kyv1,Kxh1,Kyh1,Khu1,Khv1,Khh1]=...
            uvhAssemblyIntPointImplicitSUPG(Iint,ndim,MUA,...
            bnod,hnod,unod,vnod,AGlennod,nnod,Cnod,mnod,qnod,muknod,h0nod,u0nod,v0nod,as0nod,ab0nod,as1nod,ab1nod,dadhnod,Bnod,Snod,rhonod,...
            Henod,deltanod,Hposnod,dnod,Dddhnod,...
            LSFMasknod,...
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
        
    end
    
end




%     if CtrlVar.InfoLevelCPU ;
%         UaInfo.CPUuvhAssembly=UaInfo.CPUuvhAssembly+toc(tAssembly) ;
%         UaInfo.CPUuvhAssemblyCounter=UaInfo.CPUuvhAssemblyCounter+1 ;
%     end
%

%% assemble right-hand side

Tint=sparseUA(neq,1); Fext=sparseUA(neq,1);

for Inod=1:MUA.nod
    
    
    Tint=Tint+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Tx(:,Inod),neq,1);
    Tint=Tint+sparseUA(MUA.connectivity(:,Inod)+neqx,ones(MUA.Nele,1),Ty(:,Inod),neq,1);
    Tint=Tint+sparseUA(MUA.connectivity(:,Inod)+2*neqx,ones(MUA.Nele,1),Th(:,Inod),neq,1);
    
    Fext=Fext+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Fx(:,Inod),neq,1);
    Fext=Fext+sparseUA(MUA.connectivity(:,Inod)+neqx,ones(MUA.Nele,1),Fy(:,Inod),neq,1);
    Fext=Fext+sparseUA(MUA.connectivity(:,Inod)+2*neqx,ones(MUA.Nele,1),Fh(:,Inod),neq,1);
end
%%

R=Tint-Fext;

% R=Tint-Fext;
% Tint=[Tx ; Ty ; Th] ;
% Rint=[Fx ; Fy ; Fh] ;

if ~Ronly
    
    
    
    % large memory version
    largeMemory=1;
    if largeMemory==1
        
        Iind=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1); Jind=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1);Xval=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1);
        istak=0;
        for Inod=1:MUA.nod
            for Jnod=1:MUA.nod
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Kxu(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Kxv(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Kxh(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Kyu(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Kyv(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Kyh(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Khu(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Khv(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Khh(:,Inod,Jnod);
                istak=istak+MUA.Nele;
            end
            
        end
        
        
        K=sparseUA(Iind,Jind,Xval,neq,neq);
        
        
        %%
        
    else
        Iind=zeros(9*MUA.nod*MUA.Nele,1); Jind=zeros(9*MUA.nod*MUA.Nele,1);Xval=zeros(9*MUA.nod*MUA.Nele,1);
        K=sparseUA(neq,neq);
        
        for Inod=1:MUA.nod
            istak=0;
            for Jnod=1:MUA.nod
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Kxu(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Kxv(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Kxh(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Kyu(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Kyv(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Kyh(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Khu(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Khv(:,Inod,Jnod);
                istak=istak+MUA.Nele;
                
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Khh(:,Inod,Jnod);
                istak=istak+MUA.Nele;
            end
            K=K+sparseUA(Iind,Jind,Xval,neq,neq);
        end
    end
end

%    Kxu2=Kxu; Kxv2=Kxv; Kxh2=Kxh;  Kyu2=Kyu; Kyv2=Kyv;  Kyh2=Kyh ; Khu2=Khu; Khv2=Khv ; Khh2=Khh;
%    save File2 Kxu2 Kxv2 Kxh2 Kyu2 Kyv2 Kyh2 Khu2 Khv2 Khh2

if CtrlVar.IncludeTG3uvhBoundaryTerm && CtrlVar.TG3
    %[Ktest,Rtest]=BoundaryIntegralFullyImplicitTG3(CtrlVar,MUA,h0,h,u0,v0,u,v,as0+ab0,as1+ab1,dt);
    %[K,rh]=BoundaryIntegralFullyImplicitTG3(coordinates,connectivity,Boundary,h0,h1,u0,v0,u1,v1,a0,a1,dt,CtrlVar)
    [Ktest,Rtest]=BoundaryIntegralFullyImplicitTG3(MUA.coordinates,MUA.connectivity,MUA.Boundary,F0.h,F1.h,F0.ub,F0.vb,F1.ub,F1.uv,F0.as+F0.ab,F1.as+F1.ab,CtrlVar.dt,CtrlVar);
    R=R+Rtest;
    if ~Ronly
        K=K+Ktest;
    end
end

minh=min(F1.h);


if minh<2*CtrlVar.ThickMin && CtrlVar.InfoLevelNonLinIt>1000   % if min thickness is approaching ThickMin give some information on h within NR loop
    msg=sprintf('In NRuvh loop, assembly stage: min(h) %-f \t max(h) %-g \n ',minh,max(F1.h)) ;
    fprintf(CtrlVar.fidlog,msg) ;
end

if ~Ronly
    if full(any(isnan(diag(K))))
        save TestSave  ;
        error(' NaN in K ' ) ;
    end 
end

if any(isnan(R))
    save TestSave  ;
    error("uvhMatrixAssembly:NaNinR"," NaN in R " ) ;
end
end




