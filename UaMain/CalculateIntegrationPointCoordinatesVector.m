function [xint,yint] = CalculateIntegrationPointCoordinatesVector(connectivity,coordinates,nip)
    
    %% intvec
    
    [Nele,nod]=size(connectivity); ndim=2;
    
    [points,weights]=sample('triangle',nip,ndim);
    
    coox=reshape(coordinates(connectivity,1),Nele,nod);
    cooy=reshape(coordinates(connectivity,2),Nele,nod);
    
    xint=zeros(Nele,nip) ; yint=xint;
    for Iint=1:nip                           % loop over integration points
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        xint(:,Iint)=coox*fun;
        yint(:,Iint)=cooy*fun;
        
    end
    
   
    
    
    
    