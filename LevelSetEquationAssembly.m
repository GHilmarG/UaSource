
function [UserVar,kv,rh]=LevelSetEquationAssembly(UserVar,CtrlVar,MUA,f0,c,u,v)
    
    
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
    
    
    ndim=2; dof=1; neq=dof*MUA.Nnodes;
    
    theta=CtrlVar.theta;
    dt=CtrlVar.dt; 
    tauSUPG=CalcSUPGtau(CtrlVar,MUA,u,v,dt);
    
    l=sqrt(TriAreaFE(MUA.coordinates,MUA.connectivity));
 
    speed=sqrt(u.*u+v.*v) ; V=abs(speed-c) ; 
    V=Nodes2EleMean(MUA.connectivity,V) ;
    mu=V.*l;  % looks resonable to me
    
    f0nod=reshape(f0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    unod=reshape(u(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    vnod=reshape(v(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    cnod=reshape(c(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    
    % u1nod=reshape(u1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    % v1nod=reshape(v1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    % RHS1nod=reshape(RHS1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
 
    
    % kappanod=reshape(kappa(MUA.connectivity,1),MUA.Nele,MUA.nod);
    tauSUPGnod=reshape(tauSUPG(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    
    % [points,weights]=sample('triangle',MUA.nip,ndim);
    
    
    d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
    b1=zeros(MUA.Nele,MUA.nod);
    
    % vector over all elements for each integration point
    for Iint=1:MUA.nip
        
        
        
        fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
            Deriv=MUA.Deriv(:,:,:,Iint);
            detJ=MUA.DetJ(:,Iint);
        else
            [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
        end
        
        
        
        % Deriv : Nele x dof x nod
        %  detJ : Nele
        
        % values at integration point
        
        f0int=f0nod*fun;
        
        uint=unod*fun;
        vint=vnod*fun;
        cint=cnod*fun;
        
        % u1int=u1nod*fun;
        % v1int=v1nod*fun;
        % RHSint=RHS1nod*fun;
        
        % kappaint=kappanod*fun;
        tauSUPGint=tauSUPGnod*fun;
        
        % du1dx=zeros(MUA.Nele,1); du0dx=zeros(MUA.Nele,1);
        % dv1dy=zeros(MUA.Nele,1); dv0dy=zeros(MUA.Nele,1);
        
        % derivatives at one integration point for all elements
        df0dx=zeros(MUA.Nele,1); df0dy=zeros(MUA.Nele,1);
        for Inod=1:MUA.nod
            
            
            df0dx=df0dx+Deriv(:,1,Inod).*f0nod(:,Inod);
            df0dy=df0dy+Deriv(:,2,Inod).*f0nod(:,Inod);
            
        end
        
        % Norm of gradient (NG)
        NG=sqrt(df0dx.*df0dx+df0dy.*df0dy); % at each integration point for all elements
        
        
        LT1=NG<1; 
        kappaint=1-1./NG ; %  positive diffusion 
        kappaint(LT1)=sin(2*pi*NG(LT1))./(2*pi*NG(LT1)+eps) ; % neg, then pos diffusion
        kappaint=mu.*kappaint;
        
        detJw=detJ*MUA.weights(Iint);
        
        % dt theta ( d(u1 h1)/dx    + d(v1 h1)/dy) + h1=
        %  h0+dt { (1-theta) a0+theta a1-(1-theta) (d(u0 h0)/dx+d(v0 h0)/dy}
        
        
        
        for Inod=1:MUA.nod
            
            SUPG=fun(Inod)+CtrlVar.Tracer.SUPG.Use*tauSUPGint.*(uint.*Deriv(:,1,Inod)+vint.*Deriv(:,2,Inod));
            
            SUPGdetJw=SUPG.*detJw;
            
            for Jnod=1:MUA.nod
                
                f1term=fun(Jnod).*SUPGdetJw;
                
                % fdu1dx=dt*theta*du1dx.*fun(Jnod).*SUPGdetJw;
                udf1dx=dt*theta*uint.*Deriv(:,1,Jnod).*SUPGdetJw;
                
                % fdv1dy=dt*theta*dv1dy.*fun(Jnod).*SUPGdetJw;
                vdf1dy=dt*theta*vint.*Deriv(:,2,Jnod).*SUPGdetJw;
                
                kdf1dx=dt*theta.*kappaint.*Deriv(:,1,Jnod).*Deriv(:,1,Inod).*detJw;
                kdf1dy=dt*theta.*kappaint.*Deriv(:,2,Jnod).*Deriv(:,2,Inod).*detJw;
                
                % d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+h1term+fdxu1+udxf1+fdyv1+vdyf1+kdxf1+kdyf1;
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+f1term+udf1dx+vdf1dy+kdf1dx+kdf1dy;
                
            end
            
            f0term=f0int.*SUPGdetJw;  % this is the h term, ie not \Delta h term (because the system is linear in h no need to write it in incremental form)
            
            RHSterm=-dt*cint.*NG.*SUPGdetJw;
            
            % RHS1term=dt*theta*RHS1.*SUPGdetJw;
            
            % fdxu0=-dt*(1-theta)*du0dx.*f0int.*SUPGdetJw;
            udf0dx=-dt*(1-theta)*df0dx.*uint.*SUPGdetJw;
            
            % fdyv0=-dt*(1-theta)*dv0dy.*f0int.*SUPGdetJw;
            vdf0dy=-dt*(1-theta)*df0dy.*vint.*SUPGdetJw;
            
            kappadfdx=-dt*(1-theta)*kappaint.*df0dx.*Deriv(:,1,Inod).*detJw;
            kappadfdy=-dt*(1-theta)*kappaint.*df0dy.*Deriv(:,2,Inod).*detJw;
            

            
            % b1(:,Inod)=b1(:,Inod)+f0term+RHS0term+RHS1term+fdxu0+udxf0+fdyv0+vdyf0+kappadhdx+kappadhdy;
            b1(:,Inod)=b1(:,Inod)+f0term+RHSterm+udf0dx+vdf0dy+kappadfdx+kappadfdy;
            
        end
    end
    
    % assemble right-hand side
    
    rh=sparseUA(neq,1);
    for Inod=1:MUA.nod
        rh=rh+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),b1(:,Inod),neq,1);
    end
    
    
    
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