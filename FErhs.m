function rhs=FErhs(fx,fy,coordinates,connectivity,nip)

    % calculates right-hand side in a vectorized form

    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2; dof=2; neq=dof*Nnodes;
    neqx=Nnodes ;
    
    fxnod=reshape(fx(connectivity,1),Nele,nod);   % Nele x nod
    fynod=reshape(fy(connectivity,1),Nele,nod);

    [points,weights]=sample('triangle',nip,ndim);
    Fx=zeros(Nele,nod);  Fy=zeros(Nele,nod);
    
    for Iint=1:nip
        
        
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [~,detJ]=derivVector(coordinates,connectivity,nip,Iint);
        
        fxint=fxnod*fun;
        fyint=fynod*fun;
        detJw=detJ*weights(Iint);
        
        
        for Inod=1:nod
            Fx(:,Inod)=Fx(:,Inod)+fxint.*fun(Inod).*detJw;
            Fy(:,Inod)=Fy(:,Inod)+fyint.*fun(Inod).*detJw;
        end
    end
    
    % assemble right-hand side
    
    rhs=zeros(neq,1);
    
    for Inod=1:nod
        
        rhs=rhs+sparse(connectivity(:,Inod),ones(Nele,1),Fx(:,Inod),neq,1);
        rhs=rhs+sparse(connectivity(:,Inod)+neqx,ones(Nele,1),Fy(:,Inod),neq,1);
    end
    
    
    
    
end



