function [UserVar,rh,kv,Qx,Qy]=LevelSetEquationAssemblyNR2consistentScalar(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1)

%% Level Set Equation with a source term
%
% $$ \partial \varphi/\partial t + v \cdot \nabla \varphi - \nabla \cdot (\kappa \nabla \varphi ) + c \| \nabla \varphi \| = 0$$
%
narginchk(15,53)
nargoutchk(2,4)

nOut=nargout;

ndim=2; dof=1; neq=dof*MUA.Nnodes;

theta=CtrlVar.LevelSetTheta;
dt=CtrlVar.dt;
CtrlVar.Tracer.SUPG.tau=CtrlVar.LevelSetSUPGtau;

isL=CtrlVar.LSF.L ; isP=CtrlVar.LSF.P ; isT=CtrlVar.LSF.T ; isC=CtrlVar.LSF.C; 

if  ~isfield(CtrlVar.LSF,"PG") ||  isnan(CtrlVar.LSF.PG)
    isPG=isL ;
else
    isPG=CtrlVar.LSF.PG;
end


f0nod=reshape(f0(MUA.connectivity,1),MUA.Nele,MUA.nod);
f1nod=reshape(f1(MUA.connectivity,1),MUA.Nele,MUA.nod);

u0nod=reshape(u0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
v0nod=reshape(v0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod

u1nod=reshape(u1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
v1nod=reshape(v1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod


c0nod=reshape(c0(MUA.connectivity,1),MUA.Nele,MUA.nod);
c1nod=reshape(c1(MUA.connectivity,1),MUA.Nele,MUA.nod);


if isC
    qx0nod=reshape(qx0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    qy0nod=reshape(qy0(MUA.connectivity,1),MUA.Nele,MUA.nod);

    qx1nod=reshape(qx1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    qy1nod=reshape(qy1(MUA.connectivity,1),MUA.Nele,MUA.nod);
end

d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
%b1=zeros(MUA.Nele,MUA.nod);
qx=zeros(MUA.Nele,MUA.nod);
qy=zeros(MUA.Nele,MUA.nod);
RTest=zeros(MUA.Nele,MUA.nod);

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
        
    end

    if isC  % consistent formulation 
        for Inod=1:MUA.nod
            dqx0dx=dqx1dx+Deriv(:,1,Inod).*qx0nod(:,Inod);
            dqy0dy=dqy1dy+Deriv(:,2,Inod).*qy0nod(:,Inod);

            dqx1dx=dqx1dx+Deriv(:,1,Inod).*qx1nod(:,Inod);
            dqy1dy=dqy1dy+Deriv(:,2,Inod).*qy1nod(:,Inod);
        end
    end

    
    
    
    
    % Norm of gradient (NG)
    NG0=sqrt(df0dx.*df0dx+df0dy.*df0dy+eps); % at each integration point for all elements
    NG1=sqrt(df1dx.*df1dx+df1dy.*df1dy+eps); % at each integration point for all elements
    n1x=-df1dx./NG1;  n1y=-df1dy./NG1;
    n0x=-df0dx./NG0;  n0y=-df0dy./NG0;
    
    % if gradient is very small, set normal to zero, and with it the cx and cy components
    I0=NG0< eps^2 ;
    %I1=NG1< eps^2 ;
    %n1x(I1)=0 ; n1y(I1)=0;
    n0x(I0)=0 ; n0y(I0)=0;
    
 
    cx1int=-c1int.*n1x ; cy1int=-c1int.*n1y;
    cx0int=-c0int.*n0x ; cy0int=-c0int.*n0y;
    
    tauSUPGint=CalcSUPGtau(CtrlVar,MUA.EleAreas,u0int-cx0int,v0int-cy0int,dt);
    % tauSUPGint=CalcSUPGtau(CtrlVar,MUA.EleAreas,u0int,v0int,dt);
    
    
    % I need to think about a good def for mu
    %
    % Idea :  sqrt( (u0int-cx0int).^2+(v0int-cy0int).^2)) .*sqrt(2*MUA.EleAreas) ;
    %
    
    
    switch lower(CtrlVar.LevelSetFABmu.Scale)
        
        case "constant"
            Scale=1 ;
        case "ucl"
            % Various options

            % 1)
            Scale =  sqrt( (u0int).^2+(v0int).^2) .*sqrt(2*MUA.EleAreas) ;
            % This makes the scale independent of the solution and there is no impact on the NR solution


             % 2) 
             % Scale =  sqrt( (u0int-cx0int).^2+(v0int-cy0int).^2) .*sqrt(2*MUA.EleAreas) ;
             % Here the solution will depend on f0. There is no contribution to the NR terms,  but if I solve again using backward Euler I
             % will not get the same solution if I advance both F0 and F1. This is acceptable in a transient theta=0.5 solution, but less
             % so if using a backward Euler when solving the fix point problem. 

             % 3)
             % Scale =  sqrt( (u1int-cx1int).^2+(v1int-cy1int).^2) .*sqrt(2*MUA.EleAreas) ;  % this creates a dependency of the scale on
                                                                                           % the solution which is currenlty not included in the NR terms 
            % Scale=1e5;                                                                    

        otherwise
            
            error("Ua:CaseNotFound","CtrlVar.LevelSetFABmu.Scale has an invalid value.")
    end
    
    if CtrlVar.LevelSetFABmu.Value==0 && isP && ~isT && ~isL
        CtrlVar.LevelSetFABmu.Value=1;
    end
    % This has the dimention l^2/t
    
    mu=Scale*CtrlVar.LevelSetFABmu.Value;  % This has the dimention l^2/t
 
    
    [kappaint0]=LevelSetEquationFAB(CtrlVar,NG0,mu);
    [kappaint1,dkappa]=LevelSetEquationFAB(CtrlVar,NG1,mu);
    
    
    detJw=detJ*MUA.weights(Iint);
    
    if any(~isfinite(n1x)) || any(~isfinite(kappaint1)) || any(~isfinite(dkappa))
        save TestSaveLSFnotFinite
        error("LevelSetEquationAssemblyNR2:notfinite","n0x, n1x kappa not finite")
    end
    
    
    
    
    
    nod=MUA.nod ; 
    
    for Inod=1:MUA.nod
        
        
         SUPG=CtrlVar.Tracer.SUPG.Use*tauSUPGint.*((u0int-cx0int).*Deriv(:,1,Inod)+(v0int-cy0int).*Deriv(:,2,Inod));
        % SUPG=CtrlVar.Tracer.SUPG.Use*tauSUPGint.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));
        %SUPGdetJw=SUPG.*detJw ; % if there is no advection term, set to zero, ie use Galerkin weighting
        
        if nOut>2
            for Jnod=1:nod
                
                
                TL=fun(Jnod);

                n1xDeriv1Jnod=n1x.*Deriv(:,1,Jnod);
                n1yDeriv2Jnod=n1y.*Deriv(:,2,Jnod);
                
                % (advection term)
                LL=dt*theta*(...
                    u1int.*Deriv(:,1,Jnod) + v1int.*Deriv(:,2,Jnod)...
                    -c1int.*(n1xDeriv1Jnod + n1yDeriv2Jnod)...
                    );
                
                % Pertubation term (diffusion)
                
                Plhs=dt*theta*...
                    +(kappaint1.*(Deriv(:,1,Jnod).*Deriv(:,1,Inod)+Deriv(:,2,Jnod).*Deriv(:,2,Inod)) ...
                    -NR*dkappa.*(n1xDeriv1Jnod+n1yDeriv2Jnod).*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod)));

                
                AddUp=((isT*TL+isL*LL).*(fun(Inod)+isPG.*SUPG)+isP*Plhs).*detJw ;

                
                % The dqx1dx and dqy1dy terms are calculated from the
                % previous interative solution, and therefore do not
                % depend on phi at this iteration step
                
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+AddUp;
                
            end
        end
        
        % Note, I solve: LSH  \phi  = - RHS
        
        %% Galerkin
        % (time derivative)
        TR=f1int-f0int;
        %Trhs=TR.*fun(Inod).*detJw ;
        
        % (advection+source term)
        LR= dt*(...
            +     theta*(u1int.*df1dx +v1int.*df1dy)...
            + (1-theta)*(u0int.*df0dx +v0int.*df0dy)...
            + theta*c1int.*NG1+(1-theta)*c0int.*NG0...
            );
        %Lrhs=LR.*fun(Inod).*detJw ;
        
        % Pertubation term (diffusion)
        Prhs=...
            dt*theta*kappaint1.*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod))...
            + dt*(1-theta)*kappaint0.*(df0dx.*Deriv(:,1,Inod)+df0dy.*Deriv(:,2,Inod));
        
        % Second-order form of the diffusion term
        D2=...
            - dt*  theta   * (dqx1dx+ dqy1dy)...
            - dt*(1-theta) * (dqx0dx+ dqy0dy) ;
        
%        ResidualStrongSUPGweighted=ResidualStrong.*SUPGdetJw;
        %%
        
        % qx= kappaint0.*df0dx ;
        % qu= kappaint0.*df0dy) ;
        qx(:,Inod)=qx(:,Inod)+kappaint1.*df1dx ;
        qy(:,Inod)=qy(:,Inod)+kappaint1.*df1dy ;
        
        %PG(:,Inod)=PG(:,Inod)+Prhs.*detJw;
        %LG(:,Inod)=LG(:,Inod)+Lrhs;
        %TG(:,Inod)=TG(:,Inod)+Trhs;
        %R(:,Inod)=R(:,Inod)+ResidualStrong;
 %       RSUPG(:,Inod)=RSUPG(:,Inod)+ResidualStrongSUPGweighted;
        
        RTest(:,Inod)=RTest(:,Inod)+((isT*TR+isL*LR ).*(fun(Inod)+isPG*SUPG)+isP*(Prhs+isC*D2)).*detJw; 
        
        
    end
end


Rh=sparseUA(neq,1);

for Inod=1:MUA.nod
    Rh=Rh+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),RTest(:,Inod),neq,1);
end

if isC
    Qx=sparseUA(neq,1); Qy=sparseUA(neq,1);
    for Inod=1:MUA.nod
        Qx=Qx+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),qx(:,Inod),neq,1);
        Qy=Qy+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),qy(:,Inod),neq,1);
    end
else
    Qx=[] ; Qy=[] ; 
end

% rh=isL*Lv+isP*Pv+isT*Tv+isPG*RSUPGv;
% norm(rh-Rh)/norm(rh);
rh=Rh; 


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