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
    
    theta=CtrlVar.LevelSetTheta;
   % theta=0; 
    dt=CtrlVar.dt;
    CtrlVar.Tracer.SUPG.tau=CtrlVar.LevelSetSUPGtau;
   
 
    Epsilon=CtrlVar.LevelSetEpsilon ; 

    mu=CtrlVar.LevelSetFABmu; % This had the dimention l^2/t

    
    
    switch CtrlVar.LevelSetPhase
        case "Initialisation"
            L=0 ; % The level-set equation only (i.e. without the pertubation term)
            P=1;    % P is the pertubation term
        case "Propagation"
            L=1;
            P=0;
        case "Propagation and FAB"
            L=1;
            P=1;
        otherwise
            error('safd')
    end
    

    %%
    
    f0nod=reshape(f0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    f1nod=reshape(f1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    u0nod=reshape(u0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    v0nod=reshape(v0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    
    u1nod=reshape(u1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    v1nod=reshape(v1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    
    
    c0nod=reshape(c0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    c1nod=reshape(c1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
   
    
    
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
        NG0=sqrt(df0dx.*df0dx+df0dy.*df0dy+eps); % at each integration point for all elements
        NG1=sqrt(df1dx.*df1dx+df1dy.*df1dy+eps); % at each integration point for all elements
        
        
        n1x=-df1dx./NG1;  n1y=-df1dy./NG1;
        n0x=-df0dx./NG0;  n0y=-df0dy./NG0;
        
        cx1int=-c1int.*n1x ; cy1int=-c1int.*n1y;  
        cx0int=-c0int.*n0x ; cy0int=-c0int.*n0y;
        
        tauSUPGint=CalcSUPGtau(CtrlVar,MUA,u0int-cx0int,v0int-cy1int,dt); 
        
        % I need to think about a good def for mu
        [kappaint0]=LevelSetEquationFAB(CtrlVar,NG0,mu);
        [kappaint1,dkappa]=LevelSetEquationFAB(CtrlVar,NG1,mu);
                
        detJw=detJ*weights(Iint);
        
        
        
        for Inod=1:MUA.nod
            
            
            
            SUPG=fun(Inod)+CtrlVar.Tracer.SUPG.Use*tauSUPGint.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));
            %SUPG=fun(Inod)+CtrlVar.Tracer.SUPG.Use*tauSUPGint.*((u0int-cx0int).*Deriv(:,1,Inod)+(v0int-cy0int).*Deriv(:,2,Inod));
            SUPGdetJw=SUPG.*detJw;
            
            if nargout>2
                for Jnod=1:MUA.nod
                    
                    Llhs=fun(Jnod).*SUPGdetJw...
                        +NR*dt*theta*(...
                        (u1int-cx1int).*Deriv(:,1,Jnod) + (v1int-cy1int).*Deriv(:,2,Jnod))... 
                        .*SUPGdetJw;

                     
                    Plhs=dt*theta*...
                        (kappaint1.*(Deriv(:,1,Jnod).*Deriv(:,1,Inod)+Deriv(:,2,Jnod).*Deriv(:,2,Inod)) ...
                        - NR*dkappa.*(n1x.*Deriv(:,1,Jnod)+n1y.*Deriv(:,2,Jnod)).*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod))) ...
                        .*detJw;
                    
                    
                    LHS=L*Llhs+P*Plhs; % .*SUPGdetJw ;
                    
                    Reg=Epsilon*fun(Jnod).*fun(Inod).*detJw ; 
                    
                    d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+LHS+Reg;
                    
                end
            end

            % Note, I solve: 
            %         
            %   LSH  \phi  = - RHS
            %
            Lrhs=(f1int-f0int).*SUPGdetJw...
                +    dt*theta* ((u1int-cx1int).*df1dx +(v1int-cy1int).*df1dy).*SUPGdetJw...
                + dt*(1-theta)*((u0int-cx0int).*df0dx +(v0int-cy0int).*df0dy).*SUPGdetJw; 

   % should I possibly write this term as
   %
   %    u1int.*dfdx+v1int.*df1dy + c1int.*NG1 
   %
            
            
            Prhs=...
                       dt*theta*kappaint1.*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod)).*detJw ...
                 + dt*(1-theta)*kappaint0.*(df0dx.*Deriv(:,1,Inod)+df0dy.*Deriv(:,2,Inod)).*detJw;
            
            
             Reg=Epsilon*(f1int-f0int).*fun(Inod).*detJw ; 
             
            RHS=L*Lrhs+P*Prhs+Reg ;
            
            b1(:,Inod)=b1(:,Inod)+RHS; 
            
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