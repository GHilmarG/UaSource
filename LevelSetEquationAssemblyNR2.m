function [UserVar,rh,kv]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1)
    
    
    %  df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
    %
    %
    % Treat the c norm(grad f0) term explicitly
    %
    % f1-f0  + dt theta ( v0 grad f0  )  + dt (1-theta) (v0 grad f1)  = dt  c norm(grad f0)
    % ->  f1 + dt (1-theta) v0 grad f1   = f0 + dt  c  norm(grad f0) - dt theta ( v0 grad f0  )
    %
    %
    %
    
    % norm grad f
    
    narginchk(11,11)
    
    ndim=2; dof=1; neq=dof*MUA.Nnodes;
    
    theta=CtrlVar.theta;
   % theta=0; 
    dt=CtrlVar.dt;
    CtrlVar.Tracer.SUPG.tau=CtrlVar.LevelSetSUPGtau;
    
    tauSUPG=CalcSUPGtau(CtrlVar,MUA,u0,v0,dt);  % TestIng  1
    % tauSUPG=CalcSUPGtau(CtrlVar,MUA,abs(u0-c0),v0*0 ,dt)  ;  % TestIng  2
    
    l=sqrt(MUA.EleAreas) ;
    
    speed0=sqrt(u0.*u0+v0.*v0) ; 
    % V=abs(speed0-c0) ; % TesTing 2 
    V=speed0; % TestIng 1
    V=Nodes2EleMean(MUA.connectivity,V) ; % consider doing this a integration points
    l=10e3 ; 
    mu=CtrlVar.LevelSetFAB*V.*l;  % looks resonable to me
    
    
    
    f0nod=reshape(f0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    f1nod=reshape(f1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    u0nod=reshape(u0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    v0nod=reshape(v0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    
    u1nod=reshape(u1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    v1nod=reshape(v1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    
    
    c0nod=reshape(c0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    c1nod=reshape(c1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    
    
    % kappanod=reshape(kappa(MUA.connectivity,1),MUA.Nele,MUA.nod);
    tauSUPGnod=reshape(tauSUPG(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    
    [points,weights]=sample('triangle',MUA.nip,ndim);
    
    
    d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
    b1=zeros(MUA.Nele,MUA.nod);
    
    
    if CtrlVar.LevelSetSolutionMethod=="Newton Raphson"
        NR=1;
    else
        NR=0;
    end
    
    
    % vector over all elements for each integration point
    for Iint=1:MUA.nip  % intergration points
        
        
        
        fun=shape_fun(Iint,ndim,MUA.nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
        
        
        
        
        % Deriv : Nele x dof x nod
        %  detJ : Nele
        
        % values at integration point
        
        f0int=f0nod*fun;
        f1int=f1nod*fun;
        
        u0int=u0nod*fun; v0int=v0nod*fun;
        u1int=u1nod*fun; v1int=v1nod*fun;
        c0int=c0nod*fun;
        c1int=c1nod*fun;
        
        
        tauSUPGint=tauSUPGnod*fun;
        
        
        
        % derivatives at one integration point for all elements
        df0dx=zeros(MUA.Nele,1); df0dy=zeros(MUA.Nele,1);
        df1dx=zeros(MUA.Nele,1); df1dy=zeros(MUA.Nele,1);
        for Inod=1:MUA.nod
            
            
            df0dx=df0dx+Deriv(:,1,Inod).*f0nod(:,Inod);
            df0dy=df0dy+Deriv(:,2,Inod).*f0nod(:,Inod);
            
            df1dx=df1dx+Deriv(:,1,Inod).*f1nod(:,Inod);
            df1dy=df1dy+Deriv(:,2,Inod).*f1nod(:,Inod);
        end
        
        % Norm of gradient (NG)
        NG0=sqrt(df0dx.*df0dx+df0dy.*df0dy+100*eps); % at each integration point for all elements
        NG1=sqrt(df1dx.*df1dx+df1dy.*df1dy+100*eps); % at each integration point for all elements
        
        
        n1x=-df1dx./NG1;  n1y=-df1dy./NG1;
        %n0x=-df0dx./NG0;  n0y=-df0dy./NG0;
        %c1xint=c1int.*n1x ; c1yint=c1int.*n1y;  
        %c0xint=c0int.*n0x ; c0yint=c0int.*n0y;
        
        
        [kappaint0]=LevelSetEquationFAB(CtrlVar,NG0,mu);
        [kappaint1,dkappa]=LevelSetEquationFAB(CtrlVar,NG1,mu);
                
        detJw=detJ*weights(Iint);
        
        
        
        for Inod=1:MUA.nod
            
            SUPG=fun(Inod)+CtrlVar.Tracer.SUPG.Use*tauSUPGint.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));
            SUPGdetJw=SUPG.*detJw;
            
            if nargout>2
                for Jnod=1:MUA.nod
                    
                    LHS=fun(Jnod).*SUPGdetJw...
                        +NR*dt*theta*(...
                        (u1int.*Deriv(:,1,Jnod) + v1int.*Deriv(:,2,Jnod))... 
                        +c1int.*(df1dx.*Deriv(:,1,Jnod)+df1dy.*Deriv(:,2,Jnod))./NG1)...
                        .*SUPGdetJw;

                    
                    LHSDiffusion=dt*theta*...
                        (kappaint1.*(Deriv(:,1,Jnod).*Deriv(:,1,Inod)+Deriv(:,2,Jnod).*Deriv(:,2,Inod)) ...
                        - NR*dkappa.*(n1x.*Deriv(:,1,Jnod)+n1y.*Deriv(:,2,Jnod)).*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod))) ...
                        .*detJw;
                    
                    
                    LHSterm=LHS+LHSDiffusion; % .*SUPGdetJw ;
                    
                    d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+LHSterm; % +kdf1dx+kdf1dy;
                    
                end
            end

            
            RHS=(f1int-f0int).*SUPGdetJw...
                + dt*theta*(...
                (u1int.*df1dx +v1int.*df1dy).*SUPGdetJw...
                + kappaint1.*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod)).*detJw ...
                ) ...
                + dt*(1-theta)*(...
                (u0int.*df0dx +v0int.*df0dy).*SUPGdetJw...
                + kappaint0.*(df0dx.*Deriv(:,1,Inod)+df0dy.*Deriv(:,2,Inod)).*detJw...
                ) ...
                +dt*(1-theta)*c0int.*NG0.*SUPGdetJw+dt*theta*c1int.*NG1.*SUPGdetJw...    
                ;
            
   
            
            
            b1(:,Inod)=b1(:,Inod)+RHS; % +udf0dx+vdf0dy+kappadfdx+kappadfdy;
            
        end
        
        
    end
    
    % assemble right-hand side
    
    rh=sparseUA(neq,1);
    for Inod=1:MUA.nod
        rh=rh+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),b1(:,Inod),neq,1);
    end
    
    
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