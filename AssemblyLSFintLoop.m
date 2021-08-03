function   [PG,LG,TG,R,RSUPG,qx,qy,d1d1]=AssemblyLSFintLoop(Iint,f0nod,f1nod,u0nod,u1nod,v0nod,v1nod,qx0nod,qy0nod,qx1nod,qy1nod,c0nod,c1nod,CtrlVar,NR,theta,dt,ndim,Deriv,detJ,weights,Nele,nod,points,EleAreas)
           


d1d1=zeros(Nele,nod,nod);
%b1=zeros(Nele,nod);
TG=zeros(Nele,nod);
PG=zeros(Nele,nod);
LG=zeros(Nele,nod);
R=zeros(Nele,nod);
RSUPG=zeros(Nele,nod);

qx=zeros(Nele,nod);
qy=zeros(Nele,nod);



isL=CtrlVar.LSF.L ; isP=CtrlVar.LSF.P ; isT=CtrlVar.LSF.T ;




fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points



f0int=f0nod*fun;
f1int=f1nod*fun;

u0int=u0nod*fun; v0int=v0nod*fun;
u1int=u1nod*fun; v1int=v1nod*fun;
c0int=c0nod*fun;
c1int=c1nod*fun;



% derivatives at one integration point for all elements
df0dx=zeros(Nele,1); df0dy=zeros(Nele,1);
df1dx=zeros(Nele,1); df1dy=zeros(Nele,1);
dqx0dx=zeros(Nele,1); dqy0dy=zeros(Nele,1);
dqx1dx=zeros(Nele,1); dqy1dy=zeros(Nele,1);

for Inod=1:nod
    
    
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


cx1int=-c1int.*n1x ; cy1int=-c1int.*n1y;
cx0int=-c0int.*n0x ; cy0int=-c0int.*n0y;


tauSUPGint=CalcSUPGtau(CtrlVar,EleAreas,u0int-cx0int,v0int-cy0int,dt);


switch lower(CtrlVar.LevelSetFABmu.Scale)
    
    case "constant"
        Scale=1 ;
    case "ucl"
        Scale =  sqrt( (u0int-cx0int).^2+(v0int-cy0int).^2) .*sqrt(2*EleAreas) ;
    otherwise
        
        error(adsf')
end

mu=Scale*CtrlVar.LevelSetFABmu.Value;  % This had the dimention l^2/t


[kappaint0]=LevelSetEquationFAB(CtrlVar,NG0,mu);
[kappaint1,dkappa]=LevelSetEquationFAB(CtrlVar,NG1,mu);


detJw=detJ*weights;

if any(~isfinite(n0x)) ||any(~isfinite(n1x)) || any(~isfinite(kappaint1)) || any(~isfinite(dkappa))
    save TestSaveLSFnotFinite
    error("LevelSetEquationAssemblyNR2:notfinite","n0x, n1x kappa not finite")
end


for Inod=1:nod
    
    
    SUPG=CtrlVar.Tracer.SUPG.Use*tauSUPGint.*((u0int-cx0int).*Deriv(:,1,Inod)+(v0int-cy0int).*Deriv(:,2,Inod));
    SUPGdetJw=SUPG.*detJw*isL;
    
    if nargout>2
        for Jnod=1:nod
            
            
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
                - NR*dkappa.*(n1x.*Deriv(:,1,Jnod)+n1y.*Deriv(:,2,Jnod)).*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod))) ...
                .*detJw;
            
            
            PGlhs = isT*fun(Jnod) + ...
                +isL*dt*theta*((u1int-cx1int).*Deriv(:,1,Jnod) + (v1int-cy1int).*Deriv(:,2,Jnod));
            PGlhs=SUPGdetJw.*PGlhs;
            % The dqx1dx and dqy1dy terms are calcualted from the
            % previous interative solution, and therefore do not
            % depend on phi at this iteration step
            
            d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+isL*Llhs+isP*Plhs+isT*Tlhs+PGlhs;
            
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
  
    qx(:,Inod)=qx(:,Inod)+kappaint1.*df1dx ;
    qy(:,Inod)=qy(:,Inod)+kappaint1.*df1dy ;
    
    PG(:,Inod)=PG(:,Inod)+Prhs;
    LG(:,Inod)=LG(:,Inod)+Lrhs;
    TG(:,Inod)=TG(:,Inod)+Trhs;
    R(:,Inod)=R(:,Inod)+ResidualStrong;
    RSUPG(:,Inod)=RSUPG(:,Inod)+ResidualStrongSUPGweighted;
    
    
end


