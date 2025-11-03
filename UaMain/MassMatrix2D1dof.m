function M=MassMatrix2D1dof(MUA)

% calculates the mass matrix, ie : M_{pq}=<N_p , N_q>

ndim=2;  neq=MUA.Nnodes;

%[points,weights]=sample('triangle',nip,ndim);
d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);

for Iint=1:MUA.nip
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        %        Deriv=MUA.Deriv(:,:,:,Iint);
        detJ=MUA.DetJ(:,Iint);
    else
        [~,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,MUA.points,Iint);
    end
    detJw=detJ*MUA.weights(Iint);
    
    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod
            d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+fun(Jnod).*fun(Inod).*detJw;
        end
    end
end

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1); 
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
istak=0;
for Inod=1:MUA.nod
    for Jnod=1:MUA.nod
        
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Xval(istak+1:istak+MUA.Nele)=d1d1(:,Inod,Jnod);
        istak=istak+MUA.Nele;
        
    end
end

M=sparseUA(Iind,Jind,Xval,neq,neq);
M=(M+M.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so


end



