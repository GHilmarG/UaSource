function [Dxx,Dyy]=StiffnessMatrix2D1dof(MUA)

% calculates the stiffness matrix, ie : Dxx_{pq}=<\partial_x N_p , \partial_x N_q>

neq=MUA.Nnodes;

%[points,weights]=sample('triangle',nip,ndim);
dxdx=zeros(MUA.Nele,MUA.nod,MUA.nod);
dydy=zeros(MUA.Nele,MUA.nod,MUA.nod);

for Iint=1:MUA.nip
    
    %fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    Deriv=MUA.Deriv(:,:,:,Iint); %  Nele x dof x nod
    detJ=MUA.DetJ(:,Iint);       %  Nele
    
    detJw=detJ*MUA.weights(Iint);
    
    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod
            dxdx(:,Inod,Jnod)=dxdx(:,Inod,Jnod)+Deriv(:,1,Jnod).*Deriv(:,1,Inod) .*detJw;
            dydy(:,Inod,Jnod)=dydy(:,Inod,Jnod)+Deriv(:,2,Jnod).*Deriv(:,2,Inod) .*detJw;
        end
    end
end

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1); 
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
Yval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

istak=0;
for Inod=1:MUA.nod
    for Jnod=1:MUA.nod
        
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Xval(istak+1:istak+MUA.Nele)=dxdx(:,Inod,Jnod);
        Yval(istak+1:istak+MUA.Nele)=dydy(:,Inod,Jnod);
        istak=istak+MUA.Nele;
        
    end
end

Dxx=sparseUA(Iind,Jind,Xval,neq,neq);
Dyy=sparseUA(Iind,Jind,Yval,neq,neq);
Dxx=(Dxx+Dxx.')/2 ; 
Dyy=(Dyy+Dyy.')/2 ;


end



