function [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistentNST(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1)

%% Level Set Equation without a source term
%
% $$ \partial \varphi/\partial t + (v-u) \cdot \nabla \varphi - \nabla \cdot (\kappa \nabla \varphi ) = 0 $$
%
narginchk(15,53)

ndim=2; dof=1; neq=dof*MUA.Nnodes;

theta=CtrlVar.LevelSetTheta;
dt=CtrlVar.dt;
CtrlVar.Tracer.SUPG.tau=CtrlVar.LevelSetSUPGtau;

isL=CtrlVar.LSF.L ; isP=CtrlVar.LSF.P ; isT=CtrlVar.LSF.T ;

f0nod=reshape(f0(MUA.connectivity,1),MUA.Nele,MUA.nod);
f1nod=reshape(f1(MUA.connectivity,1),MUA.Nele,MUA.nod);

u0nod=reshape(u0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
v0nod=reshape(v0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod

u1nod=reshape(u1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
v1nod=reshape(v1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod


c0nod=reshape(c0(MUA.connectivity,1),MUA.Nele,MUA.nod);
c1nod=reshape(c1(MUA.connectivity,1),MUA.Nele,MUA.nod);


qx0nod=reshape(qx0(MUA.connectivity,1),MUA.Nele,MUA.nod);
qy0nod=reshape(qy0(MUA.connectivity,1),MUA.Nele,MUA.nod);

qx1nod=reshape(qx1(MUA.connectivity,1),MUA.Nele,MUA.nod);
qy1nod=reshape(qy1(MUA.connectivity,1),MUA.Nele,MUA.nod);


d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
%b1=zeros(MUA.Nele,MUA.nod);
TG=zeros(MUA.Nele,MUA.nod);
PG=zeros(MUA.Nele,MUA.nod);
LG=zeros(MUA.Nele,MUA.nod);
R=zeros(MUA.Nele,MUA.nod);
RSUPG=zeros(MUA.Nele,MUA.nod);
qx=zeros(MUA.Nele,MUA.nod);
qy=zeros(MUA.Nele,MUA.nod);


if CtrlVar.LevelSetSolutionMethod=="Newton Raphson"
    NR=1;
else
    NR=0;
end


% vector over all elements for each  integration point
for Iint=1:MUA.nip  %Integration points
    
    
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
    
    
    f0int=f0nod*fun;
    f1int=f1nod*fun;
    
    u0int=u0nod*fun; v0int=v0nod*fun;
    u1int=u1nod*fun; v1int=v1nod*fun;
    c0int=c0nod*fun;
    c1int=c1nod*fun;
    
    
    
    % derivatives at one integration point for all elements
    df0dx=zeros(MUA.Nele,1); df0dy=zeros(MUA.Nele,1);
    df1dx=zeros(MUA.Nele,1); df1dy=zeros(MUA.Nele,1);
    dqx0dx=zeros(MUA.Nele,1); dqy0dy=zeros(MUA.Nele,1);
    dqx1dx=zeros(MUA.Nele,1); dqy1dy=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        
        df0dx=df0dx+Deriv(:,1,Inod).*f0nod(:,Inod);
        df0dy=df0dy+Deriv(:,2,Inod).*f0nod(:,Inod);
        
        df1dx=df1dx+Deriv(:,1,Inod).*f1nod(:,Inod);
        df1dy=df1dy+Deriv(:,2,Inod).*f1nod(:,Inod);
        
        dqx0dx=dqx1dx+Deriv(:,1,Inod).*qx0nod(:,Inod);
        dqy0dy=dqy1dy+Deriv(:,2,Inod).*qy0nod(:,Inod);
        
        dqx1dx=dqx1dx+Deriv(:,1,Inod).*qx1nod(:,Inod);
        dqy1dy=dqy1dy+Deriv(:,2,Inod).*qy1nod(:,Inod);
    end
    
    
    
    
    
    % Norm of gradient (NG)
    NG0=sqrt(df0dx.*df0dx+df0dy.*df0dy); % at each integration point for all elements
    NG1=sqrt(df1dx.*df1dx+df1dy.*df1dy); % at each integration point for all elements
    n1x=-df1dx./NG1;  n1y=-df1dy./NG1;
    n0x=-df0dx./NG0;  n0y=-df0dy./NG0;
    
    % if gradient is very small, set normal to zero, and with it the cx and cy components
    I0=NG0< eps^2 ;
    I1=NG1< eps^2 ;
    n1x(I1)=0 ; n1y(I1)=0;
    n0x(I0)=0 ; n0y(I0)=0;
    
    
    cx1int=c1int.*n1x ; cy1int=c1int.*n1y;
    cx0int=c0int.*n0x ; cy0int=c0int.*n0y;
    
    %% limit cx-u and cy-v where it is suffiently far away from the zero level
    
    
    
    %%
    
    tauSUPGint=CalcSUPGtau(CtrlVar,MUA.EleAreas,u0int-cx0int,v0int-cy0int,dt);
    %tauSUPGint=CalcSUPGtau(CtrlVar,MUA,u0int,v0int,dt);
    
    % I need to think about a good def for mu
    %
    % Idea :  sqrt( (u0int-cx0int).^2+(v0int-cy0int).^2)) .*sqrt(2*MUA.EleAreas) ;
    %
    
    
    switch lower(CtrlVar.LevelSetFABmu.Scale)
        
        case "constant"
            Scale=1 ;
        case "ucl"
            Scale =  sqrt( (u0int-cx0int).^2+(v0int-cy0int).^2) .*sqrt(2*MUA.EleAreas) ;
        otherwise
            
            error("Ua:CaseNotFound","CtrlVar.LevelSetFABmu.Scale has an invalid value.")
    end
    
    mu=Scale*CtrlVar.LevelSetFABmu.Value;  % This has the dimention l^2/t
    
    
    [kappaint0]=LevelSetEquationFAB(CtrlVar,NG0,mu);
    [kappaint1,dkappa]=LevelSetEquationFAB(CtrlVar,NG1,mu);
    
    
    detJw=detJ*MUA.weights(Iint);
    
    if any(~isfinite(n0x)) ||any(~isfinite(n1x)) || any(~isfinite(kappaint1)) || any(~isfinite(dkappa))
        save TestSaveLSFnotFinite
        error("LevelSetEquationAssemblyNR2:notfinite","n0x, n1x kappa not finite")
    end
    
    isPG=isL ;
    for Inod=1:MUA.nod
        
        
        SUPG=CtrlVar.Tracer.SUPG.Use*tauSUPGint.*((u0int-cx0int).*Deriv(:,1,Inod)+(v0int-cy0int).*Deriv(:,2,Inod));
        SUPGdetJw=SUPG.*detJw ; % if there is no advection term, set to zero, ie use Galerkin weighting
        
        if nargout>2
            for Jnod=1:MUA.nod
                
                
                Tlhs=fun(Jnod).*fun(Inod).*detJw;
                
                % (advection term)
                Llhs=...
                    +dt*theta*(...
                    (u1int-cx1int).*Deriv(:,1,Jnod) + (v1int-cy1int).*Deriv(:,2,Jnod))...
                    .*fun(Inod).*detJw;
                
                
                % It might appear one has forgotten  to linearize cx, but it turns out this linearisation term is equal to zero. So the only
                % contribution of the (u1int-cx1int).*df1dx term stems from the linearisation of df1dx, giving just the usual
                % u1int-cx1int).*Deriv(:,1,Jnod) + (v1int-cy1int).*Deriv(:,2,Jnod))...
                
                % Pertubation term (diffusion)
                Plhs=dt*theta*...
                    +(kappaint1.*(Deriv(:,1,Jnod).*Deriv(:,1,Inod)+Deriv(:,2,Jnod).*Deriv(:,2,Inod)) ...
                    -NR*dkappa.*(n1x.*Deriv(:,1,Jnod)+n1y.*Deriv(:,2,Jnod)).*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod))) ...
                    .*detJw;
                
                
                PGlhs = isT*fun(Jnod) + ...
                    +isL*dt*theta*((u1int-cx1int).*Deriv(:,1,Jnod) + (v1int-cy1int).*Deriv(:,2,Jnod));
                PGlhs=SUPGdetJw.*PGlhs;
                % The dqx1dx and dqy1dy terms are calculated from the
                % previous interative solution, and therefore do not
                % depend on phi at this iteration step
                
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+isL*Llhs+isP*Plhs+isT*Tlhs+isPG*PGlhs;
                
            end
        end
        
        % Note, I solve: LSH  \phi  = - RHS
        
        %% Galerkin
        % (time derivative)
        Trhs=(f1int-f0int).*fun(Inod).*detJw ;
        
        % (advection)
        Lrhs= ...
            +    dt*theta* ((u1int-cx1int).*df1dx +(v1int-cy1int).*df1dy).*fun(Inod).*detJw ...
            + dt*(1-theta)*((u0int-cx0int).*df0dx +(v0int-cy0int).*df0dy).*fun(Inod).*detJw ;
        
        % Pertubation term (diffusion)
        Prhs=...
            dt*theta*kappaint1.*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod)).*detJw ...
            + dt*(1-theta)*kappaint0.*(df0dx.*Deriv(:,1,Inod)+df0dy.*Deriv(:,2,Inod)).*detJw;
        
        
        %% Petrov
        ResidualStrong=isT*(f1int-f0int)+...
            + isL*dt*  theta   * ((u1int-cx1int).*df1dx +(v1int-cy1int).*df1dy)....
            + isL*dt*(1-theta) * ((u0int-cx0int).*df0dx +(v0int-cy0int).*df0dy)...
            - isP*dt*  theta   * (dqx1dx+ dqy1dy)...
            - isP*dt*(1-theta) * (dqx0dx+ dqy0dy) ;
        
        
        ResidualStrongSUPGweighted=ResidualStrong.*SUPGdetJw;
        %%
        
        % qx= kappaint0.*df0dx ;
        % qu= kappaint0.*df0dy) ;
        qx(:,Inod)=qx(:,Inod)+kappaint1.*df1dx ;
        qy(:,Inod)=qy(:,Inod)+kappaint1.*df1dy ;
        
        PG(:,Inod)=PG(:,Inod)+Prhs;
        LG(:,Inod)=LG(:,Inod)+Lrhs;
        TG(:,Inod)=TG(:,Inod)+Trhs;
        R(:,Inod)=R(:,Inod)+ResidualStrong;
        RSUPG(:,Inod)=RSUPG(:,Inod)+ResidualStrongSUPGweighted;
        
        
    end
end

Pv=sparseUA(neq,1);
Lv=sparseUA(neq,1);
Tv=sparseUA(neq,1);
Qx=sparseUA(neq,1);
Qy=sparseUA(neq,1);
Rv=sparseUA(neq,1);
RSUPGv=sparseUA(neq,1);

for Inod=1:MUA.nod
    Pv=Pv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),PG(:,Inod),neq,1);
    Lv=Lv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),LG(:,Inod),neq,1);
    Tv=Tv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),TG(:,Inod),neq,1);
    Rv=Rv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),R(:,Inod),neq,1);
    RSUPGv=RSUPGv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),RSUPG(:,Inod),neq,1);
    Qx=Qx+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),qx(:,Inod),neq,1);
    Qy=Qy+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),qy(:,Inod),neq,1);
end

rh=isL*Lv+isP*Pv+isT*Tv+isPG*RSUPGv;

if nargout>2
    Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1); Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
    istak=0;
    
    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod
            Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
            Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
            Xval(istak+1:istak+MUA.Nele)=d1d1(:,Inod,Jnod);
            istak=istak+MUA.Nele;
        end
    end
    
    kv=sparseUA(Iind,Jind,Xval,neq,neq);
end

end