function P=NodalFormFunctionInfluence(MUA)


% calculates:  M_p}=<N_p , 1>

ndim=2;  neq=MUA.Nnodes;

%[points,weights]=sample('triangle',nip,ndim);
d1d1=zeros(MUA.Nele,MUA.nod);

for Iint=1:MUA.nip
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        %        Deriv=MUA.Deriv(:,:,:,Iint);
        detJ=MUA.DetJ(:,Iint);
    else
        [~,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    detJw=detJ*MUA.weights(Iint);
    
    for Inod=1:MUA.nod
        %for Jnod=1:MUA.nod
            d1d1(:,Inod)=d1d1(:,Inod)+fun(Inod).*detJw;
            %d1d1(:,Inod)=d1d1(:,Inod)+detJw;
        %end
    end
end

Iind=zeros(MUA.nod*MUA.Nele,1); 
%Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
Xval=zeros(MUA.nod*MUA.Nele,1);
istak=0;
for Inod=1:MUA.nod
    %for Jnod=1:MUA.nod
        
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
     %   Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Xval(istak+1:istak+MUA.Nele)=d1d1(:,Inod);
        istak=istak+MUA.Nele;
        
    %end
end

Jind=Iind*0+1;
P=sparseUA(Iind,Jind,Xval,neq,1);
P=full(P);
%M=(M+M.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so


end
