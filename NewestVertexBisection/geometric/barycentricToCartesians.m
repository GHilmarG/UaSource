function cartesians = barycentricToCartesians(mesh, barycentric)
%BARYCENTRICTOCARTESIANS    Single barycentric coordinate to multiple
%cartesian coordinates
%   CARTESIANS = BARYCENTRICTOCARTESIANS(MESH, BARYCENTRIC) yields a
%   nE-by-dimSpace array containing the points obtained using the given
%   barycentric coordinate as local coordinate on all the elements.
%   Works for arbitrary-dimensional meshes.
%
%   Author: Josef Kemetmueller - 16.12.2013

assert(length(barycentric)==numel(barycentric), ...
                 'Barycentric coordinate must be given as a single vector');
nE = numElements(mesh);
cartesians = zeros(nE, dimSpace(mesh));
for d = 1:dimMesh(mesh)+1
    cartesians = cartesians + barycentric(d)*mesh.coordinates(mesh.elements(:,d),:);
end