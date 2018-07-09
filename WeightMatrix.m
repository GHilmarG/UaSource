function M=WeightMatrix(coordinates,connectivity,nip)
    
   
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	ndim=2;  neq=Nnodes;
	
    [points,weights]=sample('triangle',nip,ndim);
    
    for Iint=1:nip
        
        
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [~,detJ]=derivVector(coordinates,connectivity,nip,Iint);
        % Deriv : Nele x dof x nod
        %  detJ : Nele
        
        
        detJw=detJ*weights(Iint);
        
        Sele=zeros(Nele,nod,nod);
        
        for Inod=1:nod
            for Jnod=1:nod
                
                Sele(:,Inod,Jnod)=Sele(:,Inod,Jnod)+fun(Jnod).*fun(Inod).*detJw ;
                
            end
        end
    end
    
    
    
    Iind=zeros(nod*nod*Nele,1); Jind=zeros(nod*nod*Nele,1); Xval=zeros(nod*nod*Nele,1);
    istak=0;
    for Inod=1:nod
        for Jnod=1:nod
            
            Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=Sele(:,Inod,Jnod);
            istak=istak+Nele;
           
        end
        
    end
    
    %%
    
    M=sparse(Iind,Jind,Xval,neq,neq);
    
    
end



