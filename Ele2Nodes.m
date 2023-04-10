function [E2N,ElePerNode] = Ele2Nodes(connectivity,Nnodes)

%%
% [E2N,ElePerNode] = Ele2Nodes(connectivity,Nnodes)
% Creates a sparse matrix, E2N, that can be used to map element values defined for each element to nodes
%
%   f = E2N * g
%
% f : nodal values calculated by averaging over elements connected to that node
%
% g :  a vector with one value per element
%
% vectorized
%
%
% See also; Nodes2EleMean
%%


[Nele,nod]=size(connectivity);

% R(p,q)=1 for connectivity(p,q)
% R is Nnodes x Nele

I= connectivity(:);
J=repmat((1:Nele)',nod,1);
E2N=sparse(I,J,1,Nnodes,Nele);

ElePerNode=sum(E2N,2);  % this is the number of elements that each node is attaced to


% divide each line by the number of elements to which corresponding node belongs to


T=spdiags(1./ElePerNode,0,Nnodes,Nnodes)*E2N;
ind=find(E2N) ;
E2N(ind)=E2N(ind).*T(ind);


ElePerNode=full(ElePerNode);

%%

end
