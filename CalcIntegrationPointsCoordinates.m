function [xint,yint] = CalcIntegrationPointsCoordinates(coordinates,connectivity,nip)
    
    % Calculates coordinates of element integration points in a vectorized way
    % (The preferred m file to use)
    %  Returns a vector (OOps changed this on 31 July, now returns an array 
   
    [Nele,nod]=size(connectivity); ndim=2;
    [points]=sample('triangle',nip,ndim);
   
    coox=reshape(coordinates(connectivity,1),Nele,nod);
    cooy=reshape(coordinates(connectivity,2),Nele,nod);
    
    %xint=zeros(nip,Nele) ; yint=zeros(nip,Nele);       
	xint=zeros(Nele,nip) ; yint=zeros(Nele,nip);     % changed on July 31  
    for Iint=1:nip
         fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
         xint(:,Iint)=coox*fun;
         yint(:,Iint)=cooy*fun;
         
    end
   % xint=xint(:); yint=yint(:);  % changed on July 31  
    
end