function [xint,yint] = CalculateIntegrationPointCoordinates(connectivity,coordinates,nip)
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    
    funInt=cell(1,3); derInt=cell(1,3); ndim=2;
    [points,weights]=sample('triangle',nip,ndim);
    
    for Iint=1:nip
        funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
    end
    
    
    xint=zeros(Nele,nip) ; yint=xint;
    for Iele=1:Nele                           % loop over elements
        
        con=connectivity(Iele,:);  % nodes of element
        coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
       
        for Iint=1:nip                           % loop over integration points
            
            
            fun=funInt{Iint} ;
            xint(Iele,Iint)=coo(:,1)'*fun ;
            yint(Iele,Iint)=coo(:,2)'*fun ;
            
        end
    end
    
    
end

