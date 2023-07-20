function P=InnerProduct_FormFunctions_with_EleIntegrationPointVariable(MUA,Fint)

% calculates T_q=<Fint,n_q> where Fint is defined at integration points and n_q are the form functions.
% Fint must have the dimensions Nele x nip
% where Nele is the number of elements and nip the number of integration points
%        
% On output P has the dimentions MUA.Nnodes x 1
%

[n1,n2]=size(Fint);

if n1~=MUA.Nele || n2~=MUA.nip
    fprintf('Fint must have the dimensions %i x %i but has on input the dimensions %i x %i \n',MUA.Nele,MUA.nip,n1,n2)
    P=[];
    return
end

ndim=2;
P=sparseUA(MUA.Nnodes,1);
R=zeros(MUA.Nele,MUA.nod);

for Iint=1:MUA.nip
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
%        Deriv=MUA.Deriv(:,:,:,Iint);
        detJ=MUA.DetJ(:,Iint);
    else
        [~,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,MUA.points,Iint);
    end
    detJw=detJ*MUA.weights(Iint);
    
    %        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    %       [Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
    
    % Deriv : Nele x dof x nod
    %  detJ : Nele
    
    % values at integration this point

    for Inod=1:MUA.nod
        R(:,Inod)=R(:,Inod)+Fint(:,Iint).*fun(Inod).*detJw;
    end
end

for Inod=1:MUA.nod
    P=P+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),R(:,Inod),MUA.Nnodes,1);
end

end

