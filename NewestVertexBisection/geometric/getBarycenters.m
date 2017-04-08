function CC = getBarycenters(mesh)
%GETBARYCENTERS    Computes barycenters of simplices of given mesh.
%   CENTERS = GETBARYCENTERS(MESH) yields the barycenters of the simplices
%   defined by mesh.elements.
%   Works for arbitrary-dimensional meshes.
%
%   Author: Josef Kemetmueller - 20.03.2013

nE = numElements(mesh);

CC = zeros(nE, dimSpace(mesh));

for n = 1:dimMesh(mesh)+1
    CC = CC+mesh.coordinates(mesh.elements(:,n),:)/(dimMesh(mesh)+1);
end

end