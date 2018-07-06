function [dIdC]=dIdCqTest(u,v,lx,ly,S,B,h,connectivity,coordinates,nip,C,m,rho,rhow,CtrlVar)
    
    
    % Calculates the gradient of the cost function with respect to C
    % This done by first calculating the gradient at nodal points, and then
    % evaluationg the resulting nodal field at integration points
    
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2;
    
    hf=(S-B)*rhow/rho;
    kH=CtrlVar.kH;
    He = HeavisideApprox(kH,h-hf);
     
    
    
    Cgrad= (1/m)*He.*(C+CtrlVar.CAdjointZero).^(-1/m-1).*(sqrt(u.*u+v.*v+CtrlVar.SpeedZero^2)).^(1/m-1).*(u.*lx+v.*ly);
    
    Cgradnod=reshape(Cgrad(connectivity,1),Nele,nod);
    
    [points,weights]=sample('triangle',nip,ndim);
    
    
    T=zeros(Nele,nod);
    
    for Iint=1:nip
        
        
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [~,detJ]=derivVector(coordinates,connectivity,nip,Iint);
        
        
        Cgradint=Cgradnod*fun;
        detJw=detJ*weights(Iint);
        
        
        for Inod=1:nod
            
            %T(:,Inod)=T(:,Inod)+Cgradint.*(uint.*lxint+vint.*lyint).*fun(Inod).*detJw;
            T(:,Inod)=T(:,Inod)+Cgradint.*fun(Inod).*detJw;
            
        end
    end
    
    dIdC=zeros(Nnodes,1);
    
    for Inod=1:nod
        dIdC=dIdC+sparse(connectivity(:,Inod),ones(Nele,1),T(:,Inod),Nnodes,1);
    end
    
end



