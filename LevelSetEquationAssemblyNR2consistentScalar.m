function [UserVar,Rh,K,Qx,Qy]=LevelSetEquationAssemblyNR2consistentScalar(UserVar,CtrlVar,MUA,F0,F1) 



f0=F0.LSF ; c0=F0.c ; u0=F0.ub ; v0=F0.vb ; f1=F1.LSF ; c1=F1.c ; u1=F1.ub ; v1=F1.vb ; qx0=F0.LSFqx ; qy0=F0.LSFqy ; qx1=F1.LSFqx ; qy1=F1.LSFqy;


% [UserVar,f0,K,Qx,Qy]=LevelSetEquationAssemblyNR2consistentScalar(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1)

%% Level Set Equation with a source term
%
% $$ \partial \varphi/\partial t + v \cdot \nabla \varphi - \nabla \cdot (\kappa \nabla \varphi ) + c \| \nabla \varphi \| = 0$$
%
%
%     K d\varphi = - f0
%
%


narginchk(5,5)
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



f1nod=reshape(f1(MUA.connectivity,1),MUA.Nele,MUA.nod);  % this is the level set function
f0nod=reshape(f0(MUA.connectivity,1),MUA.Nele,MUA.nod);  % this is the level set function

u0nod=reshape(u0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
v0nod=reshape(v0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod

u1nod=reshape(u1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
v1nod=reshape(v1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod

coox=reshape(MUA.coordinates(MUA.connectivity,1),MUA.Nele,MUA.nod);
cooy=reshape(MUA.coordinates(MUA.connectivity,2),MUA.Nele,MUA.nod);


if  ~contains(CtrlVar.CalvingLaw.Evaluation,"-int-")
    c0nod=reshape(c0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    c1nod=reshape(c1(MUA.connectivity,1),MUA.Nele,MUA.nod);
end

% additonal variables for sliding law evaluation at int point

h1nod=reshape(F1.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
h0nod=reshape(F0.h(MUA.connectivity,1),MUA.Nele,MUA.nod);

s1nod=reshape(F1.s(MUA.connectivity,1),MUA.Nele,MUA.nod);
s0nod=reshape(F0.s(MUA.connectivity,1),MUA.Nele,MUA.nod);

rhonod=reshape(F0.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);



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
    
    xint=coox*fun;  % coordinates of this integration point for all elements
    yint=cooy*fun;


    f0int=f0nod*fun;
    f1int=f1nod*fun;
    
    u0int=u0nod*fun; v0int=v0nod*fun;
    u1int=u1nod*fun; v1int=v1nod*fun;
  
    % additonal variables for sliding law evaluation at int point
    h1int=h1nod*fun;
    h0int=h0nod*fun;
    s1int=s1nod*fun;
    s0int=s0nod*fun;
    rhoint=rhonod*fun;
    
    % derivatives at one integration point for all elements
    df0dx=zeros(MUA.Nele,1); df0dy=zeros(MUA.Nele,1);
    df1dx=zeros(MUA.Nele,1); df1dy=zeros(MUA.Nele,1);
    dqx0dx=zeros(MUA.Nele,1); dqy0dy=zeros(MUA.Nele,1);
    dqx1dx=zeros(MUA.Nele,1); dqy1dy=zeros(MUA.Nele,1);

    exx1=zeros(MUA.Nele,1);
    eyy1=zeros(MUA.Nele,1);
    exy1=zeros(MUA.Nele,1);
    exx0=zeros(MUA.Nele,1);
    eyy0=zeros(MUA.Nele,1);
    exy0=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod


        df0dx=df0dx+Deriv(:,1,Inod).*f0nod(:,Inod);
        df0dy=df0dy+Deriv(:,2,Inod).*f0nod(:,Inod);

        df1dx=df1dx+Deriv(:,1,Inod).*f1nod(:,Inod);
        df1dy=df1dy+Deriv(:,2,Inod).*f1nod(:,Inod);

        exx1=exx1+Deriv(:,1,Inod).*u1nod(:,Inod);
        eyy1=eyy1+Deriv(:,2,Inod).*v1nod(:,Inod);
        exy1=exy1+0.5*(Deriv(:,1,Inod).*v1nod(:,Inod) + Deriv(:,2,Inod).*u1nod(:,Inod));
       

        exx0=exx0+Deriv(:,1,Inod).*u0nod(:,Inod);
        eyy0=eyy0+Deriv(:,2,Inod).*v0nod(:,Inod);
        exy0=exy0+0.5*(Deriv(:,1,Inod).*v0nod(:,Inod) + Deriv(:,2,Inod).*u0nod(:,Inod));


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
 
    n0x(I0)=0 ; n0y(I0)=0;




    if  contains(CtrlVar.CalvingLaw.Evaluation,"-int-")



        nArgs=nargin('DefineCalvingAtIntegrationPoints');

        if nArgs==5

            F1int.ub=u1int;
            F1int.vb=v1int;
            F1int.h=h1int;
            F1int.s=s1int;
            F1int.rho=rhoint;
            F1int.exx=exx1;
            F1int.eyy=eyy1;
            F1int.exy=exy1;
            F1int.x=xint;
            F1int.y=yint;
            F1int.S=F1.S(1);
            F1int.rhow=F1.rhow;

            F0int.ub=u0int;
            F0int.vb=v0int;
            F0int.h=h0int;
            F0int.s=s0int;
            F0int.rho=rhoint;
            F0int.exx=exx0;
            F0int.eyy=eyy0;
            F0int.exy=exy0;
            F0int.x=xint;
            F0int.y=yint;
            F0int.S=F0.S(1);
            F0int.rhow=F0.rhow;


            [c1int,dcDdfdx1,dcDdfdy1]=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,df1dx,df1dy,F1int) ;
            [c0int]=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,df0dx,df0dy,F0int) ;

        else

            [c1int,dcDdfdx1,dcDdfdy1]=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,df1dx,df1dy,u1int,v1int,h1int,s1int,F1.S(1),xint,yint) ;
            c0int=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,df0dx,df0dy,u0int,v0int,h0int,s0int,F0.S(1),xint,yint) ;

        end

    else
        c0int=c0nod*fun;
        c1int=c1nod*fun;
        dcDdfdx1=0;
        dcDdfdy1=0;
    end

   % cx1int=c1int.*n1x ; cy1int=c1int.*n1y;
    cx0int=c0int.*n0x ; cy0int=c0int.*n0y;  % only used when calculating tauSUPG

    % N=20 ; [u0int(1:N).*n0x(1:N)+v0int(1:N).*n0y(1:N) c0int(1:N)]
    %
    % If calving rate is set equal to ice-velocity normal to calving front then these should be zero:
    %
    % norm(u0int.*n0x+v0int.*n0y-c0int)  
    % norm(u0int.*n0x+v0int.*n0y-(cx0int.*n0x+cy0int.*n0y) )
    %
    % N=4 ; [u1int(1:N).*df1dx(1:N)+v1int(1:N).*df1dy(1:N) c1int(1:N).*NG1(1:N) (cx1int(1:N).*n1x(1:N)+cy1int(1:N).*n1y(1:N)).*NG1(1:N)]
    %%


     


    % let 
    %
    %   c0int.*n0x+v0int.*n0y-c0int
    %
    % be the speed scale. This is zero in the limit when calving rate equals veocity normal to the front
    % where the advective terms becomes zero.
    %
    % Using 
    %
    %   (u0int-cx0int,v0int-cy0int)
    %
    % does in this case not equal exactly zero.
    %
    % 

    % tauSUPGint=CalcSUPGtau(CtrlVar,MUA.EleAreas,u0int-cx0int,v0int-cy0int,dt);

    % CtrlVar.Tracer.SUPG.tau='tau2' ;   %  inversly weighted average of spatial and temporal tau. This is dt/2 if speed -> 0
    CtrlVar.Tracer.SUPG.tau='tau1' ;   %   typical textbook recomendation for spatially constant (and non-zero) speed for linear advection equation. This is dt/6 if speed -> 0
    tauSUPGint=CalcSUPGtau(CtrlVar,MUA.EleAreas,u0int.*n0x+v0int.*n0y-c0int,0,dt); 
   


    if CtrlVar.LevelSetFABmu.Value==0 && isP && ~isT && ~isL
        CtrlVar.LevelSetFABmu.Value=1;
    end

    % This has the dimention l^2/t
    % I need to think about a good def for mu
    %
    % Idea :  sqrt( (u0int-cx0int).^2+(v0int-cy0int).^2)) .*sqrt(2*MUA.EleAreas) ;
    %

    % The Scale (or \mu in the UaCompendium) has the units m^2/t
    %
    % Various options :
    switch lower(CtrlVar.LevelSetFABmu.Scale)

        case "constant"

            %  When solvingt the fix-point problem, I just want the diffusion term to be driven to zero This will happen once the min of
            %  the potential is reached, and the minimum is independent of mu. So, arguably, one should just set mu to some fixed
            %  constant value to ensure that ||\grad \varphi||= 0 at every point within the domain.

            Scale =1 ;
            mu=Scale*CtrlVar.LevelSetFABmu.Value;

        case {"u-cl","-u-cl-"}

            % 1)
            Scale =  sqrt( (u0int-cx0int).^2+(v0int-cy0int).^2) .*sqrt(2*MUA.EleAreas) ;
            % Here the solution will depend on f0. There is no contribution to the NR terms,  but if I solve again using backward Euler I
            % will not get the same solution if I advance both F0 and F1. This is acceptable in a transient theta=0.5 solution, but less
            % so if using a backward Euler when solving the fix point problem.
            mu=Scale*CtrlVar.LevelSetFABmu.Value;
        case {"ul","ucl","-ucl-","-ul-"}   % I should not include ucl here, but this is to make it compatiable with earlier inputs


            % 2)
            Scale =  sqrt( (u0int).^2+(v0int).^2) .*sqrt(2*MUA.EleAreas) ;
            % This makes the scale independent of the solution and there is no impact on the NR solution
            % but in the particular case where v=c  the diffusion terms does not drop out. This might potentially cause a shift in the
            % steady-state solution for (xc,yc)

            % 3)
            % Scale =  sqrt( (u1int-cx1int).^2+(v1int-cy1int).^2) .*sqrt(2*MUA.EleAreas) ;  % this creates a dependency of the scale on
            % the solution which is currenlty not included in the NR terms
            % Scale=1e5;
            mu=Scale*CtrlVar.LevelSetFABmu.Value;
        otherwise

            error("Ua:CaseNotFound","CtrlVar.LevelSetFABmu.Scale has an invalid value.")
    end


    
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
               
        if nOut>2
            for Jnod=1:nod
                
                
                TL=fun(Jnod);

                n1xDeriv1Jnod=n1x.*Deriv(:,1,Jnod);
                n1yDeriv2Jnod=n1y.*Deriv(:,2,Jnod);

                LL=dt*theta*(...
                    u1int.*Deriv(:,1,Jnod) + v1int.*Deriv(:,2,Jnod)...
                    +c1int.*(df1dx.*Deriv(:,1,Jnod)+df1dy.*Deriv(:,2,Jnod))./(sqrt(df1dx.*df1dx+df1dy.*df1dy+eps))...
                    +(dcDdfdx1.*Deriv(:,1,Jnod)+dcDdfdy1.*Deriv(:,2,Jnod)).*NG1 ...
                    );
                

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
        
        
        % (advection+source term)
        LR= dt*(...
            +     theta*(u1int.*df1dx +v1int.*df1dy)...
            + (1-theta)*(u0int.*df0dx +v0int.*df0dy)...
            + theta*c1int.*NG1+(1-theta)*c0int.*NG0...      % NG1=sqrt(df1dx.*df1dx+df1dy.*df1dy)
            );
 
        
        % Pertubation term (diffusion)
        Prhs=...
            dt*theta*kappaint1.*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod))...
            + dt*(1-theta)*kappaint0.*(df0dx.*Deriv(:,1,Inod)+df0dy.*Deriv(:,2,Inod));
        
        % Second-order form of the diffusion term
        D2=...
            - dt*  theta   * (dqx1dx+ dqy1dy)...
            - dt*(1-theta) * (dqx0dx+ dqy0dy) ;
        

        qx(:,Inod)=qx(:,Inod)+kappaint1.*df1dx ;
        qy(:,Inod)=qy(:,Inod)+kappaint1.*df1dy ;
        
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
% f0=Rh; 


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
    
    K=sparseUA(Iind,Jind,Xval,neq,neq);
end

if ~isL
   
    % if I'm not including the advection term (ie only the transient and/or the diffusion term)
    % the system is symmetric

    K=(K+K')/2 ;

end
   









end