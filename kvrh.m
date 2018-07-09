
function [kv,rh]=kvrh(s,h,coordinates,connectivity,nip,etaInt,gfint,beta2,alpha,rho,rhow,g)
    
    % calculates the FE matrix and right-hand side in a vectorized form
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2; dof=2; neq=dof*Nnodes;
    neqx=Nnodes ; 
    
    hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
    snod=reshape(s(connectivity,1),Nele,nod);
    %unod=reshape(u(connectivity,1),Nele,nod);
    %vnod=reshape(v(connectivity,1),Nele,nod);
    
    beta2nod=reshape(beta2(connectivity,1),Nele,nod);
    
    
    rhog=rho*g; ca=cos(alpha); sa=sin(alpha);
   
    
    [points,weights]=sample('triangle',nip,ndim);
    % get local coordinates and weights
    
    
    kv=sparse(neq,neq); 
    d1d1=zeros(Nele,nod,nod); d2d2=zeros(Nele,nod,nod);  d1d2=zeros(Nele,nod,nod);  
    b1=zeros(Nele,nod);  b2=zeros(Nele,nod);
    
    for Iint=1:nip
        
        
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
        % Deriv : Nele x dof x nod
        %  detJ : Nele
        
        % values at integration this point
        hint=hnod*fun;
        %sint=snod*fun;
        %uint=unod*fun;
        %vint=vnod*fun;
        beta2int=beta2nod*fun; beta2int=beta2int.*gfint(:,Iint);
        etaint=etaInt(:,Iint) ;  % I could consider calculating this here
        
        
        % derivatives at this integration point for all elements
        dsdx=zeros(Nele,1); dhdx=zeros(Nele,1);
        dsdy=zeros(Nele,1); dhdy=zeros(Nele,1);
        for Inod=1:nod
            
            dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
            dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
            dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
            dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
            
        end
        
        detJw=detJ*weights(Iint);
        
        
        for Inod=1:nod
            for Jnod=1:nod
                
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+etaint.*hint.*(4*Deriv(:,1,Inod).*Deriv(:,1,Jnod)+Deriv(:,2,Inod).*Deriv(:,2,Jnod)).*detJw;
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+beta2int.*fun(Jnod).*fun(Inod).*detJw;
                
                d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)+etaint.*hint.*(4*Deriv(:,2,Inod).*Deriv(:,2,Jnod)+Deriv(:,1,Inod).*Deriv(:,1,Jnod)).*detJw;
                d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)+beta2int.*fun(Jnod).*fun(Inod).*detJw;
                
                d1d2(:,Inod,Jnod)=d1d2(:,Inod,Jnod)+etaint.*hint.*(2*Deriv(:,1,Inod).*Deriv(:,2,Jnod)+Deriv(:,2,Inod).*Deriv(:,1,Jnod)).*detJw;
                %d2d1(:,Inod,Jnod)=d2d1(:,Inod,Jnod)+etaint.*hint.*(2*Deriv(:,2,Inod).*Deriv(:,1,Jnod)+Deriv(:,1,Inod).*Deriv(:,2,Jnod)).*detJw;
                
            end
            
            t1=rhog*hint.*((dsdx-(1-rho/rhow).*dhdx)*ca-sa).*fun(Inod);
            %t2=-0.5*rhog.*(1-rho/rhow).*hint.^2.*Deriv(:,1,Inod); $ wrong
            t2=-0.5*ca*rhog.*(1-rho/rhow).*hint.^2.*Deriv(:,1,Inod); % correction on 11 march 2010
            b1(:,Inod)=b1(:,Inod)-(t1+t2).*detJw;
            
            t1=rhog.*hint.*(dsdy-(1-rho/rhow).*dhdy).*ca.*fun(Inod);
            %t2=-0.5*rhog.*(1-rho/rhow).*hint.^2.*Deriv(:,2,Inod);
            t2=-0.5*ca*rhog.*(1-rho/rhow).*hint.^2.*Deriv(:,2,Inod); % correction on 11 march 2010
            b2(:,Inod)=b2(:,Inod)-(t1+t2).*detJw;
            
        end
    end
    
    d2d1=d1d2;  % because I know that this term is symmetric
    
    
    % assemble right-hand side
    
    rh=sparse(neq,1);
    for Inod=1:nod
        
        rh=rh+sparse(connectivity(:,Inod),ones(Nele,1),b1(:,Inod),neq,1);
        rh=rh+sparse(connectivity(:,Inod)+neqx,ones(Nele,1),b2(:,Inod),neq,1);
        
    end
    
    rh=full(rh);
    
    % assemble matrix
    iSparse=1;  % actually there seems to be no difference between these two approaches
    if iSparse==1
        % uses the sparse function less often
        
        
        Iind=zeros(nod*Nele*4,1); Jind=zeros(nod*Nele*4,1);Xval=zeros(nod*Nele*4,1);
        for Inod=1:nod
            istak=0;
            for Jnod=1:nod
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=d1d1(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=d2d2(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=d1d2(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=d2d1(:,Inod,Jnod);
                istak=istak+Nele;
                
            end
            kv=kv+sparse(Iind,Jind,Xval,neq,neq);
            
            
        end

    else
        % creates the sparse matrix in steps, requires no extra
        for Inod=1:nod
            for Jnod=1:nod
                kv=kv+sparse(connectivity(:,Inod),connectivity(:,Jnod),d1d1(:,Inod,Jnod),neq,neq);
                kv=kv+sparse(connectivity(:,Inod)+neqx,connectivity(:,Jnod)+neqx,d2d2(:,Inod,Jnod),neq,neq);
                kv=kv+sparse(connectivity(:,Inod),connectivity(:,Jnod)+neqx,d1d2(:,Inod,Jnod),neq,neq);
                kv=kv+sparse(connectivity(:,Inod)+neqx,connectivity(:,Jnod),d2d1(:,Inod,Jnod),neq,neq);
                
            end
        end
        
    end
    
 
    
    kv=(kv+kv.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
                     % Note: for numerical verificatin of distributed parameter gradient it is important to
                     % not to use the complex conjugate transpose.
    
end



