 
function [J11,J12,J21,J22,detJ,Deriv]=deriv(coordinates,connectivity,nip,Iint)

    [Nele,nod]=size(connectivity); ndim=2; dof=2;
    [points,weights]=sample('triangle',nip,ndim);
    J11=zeros(Nele,1) ; J12=J11 ; J21=J11 ; J22=J11;
    iJ11=zeros(Nele,1) ; iJ12=iJ11 ; iJ21=iJ11 ; iJ22=iJ11;
    
    detJ=J11; 
    
    Deriv=zeros(Nele,dof,nod);
    
for Iele=1:Nele
    % gather local quantities from global arrays
    % note the nodal numbering is clockwise!
    con=connectivity(Iele,:);  % nodes of element
    coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
    
    
    %for Iint=1:nip                           % loop over integration points
        
        der=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
        J=der*coo;
        J11(Iele,Iint)=J(1,1);
        J12(Iele,Iint)=J(1,2);
        J21(Iele,Iint)=J(2,1);
        J22(Iele,Iint)=J(2,2);
        iJ=inv(J);
        iJ11(Iele,Iint)=iJ(1,1);
        iJ12(Iele,Iint)=iJ(1,2);
        iJ21(Iele,Iint)=iJ(2,1);
        iJ22(Iele,Iint)=iJ(2,2);
        
        detJ(Iele,Iint)=det(J);  % det(dof x dof) matrix
        Deriv(Iele,:,:)=J\der;
        
   % end
    
    
    
end


