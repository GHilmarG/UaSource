function element2neighbors = getNeighbors(mesh)
%GETNEIGHBORS    Returns the indices of the neighboring elements
%   ELEMENT2NEIGHBORS = GETNEIGHBORS(MESH) returns a nE-by-(dimMesh+1)
%   array containing the row indices of the neighbors of the elements.
%   ELEMENT2NEIGHBORS(i,j) is the neighbor of the i-th element
%   corresponding to its j-th face: elements(i,(1:dimMesh+1)~=j).
%
%   Works for arbitrary-dimensional meshes. 
%
% Author: Josef Kemetmueller - 20.03.14

nE = numElements(mesh);
% Generate the hyperfaces
% Get hyperfaceorder-matrix: [~1],[~2],[~3],[~4],...
hyperfaceorder = simplexHyperfaces(dimMesh(mesh));
% Generate the hyperfaces
hyperfaces = sort(reshape(mesh.elements(:,hyperfaceorder),[],dimMesh(mesh)),2);
% Find the corresponding faces
[~,I,J] = unique(hyperfaces,'rows');
A = I(J);
half1 = find(A~=(1:(dimMesh(mesh)+1)*nE)');
half2 = A(half1);
element2neighbors = zeros(nE,dimMesh(mesh)+1);
if (nE~=1)
    element2neighbors(half2) = rem(half1-1,nE)+1;
    element2neighbors(half1) = rem(half2-1,nE)+1;
end