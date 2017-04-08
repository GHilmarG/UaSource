function [mesh, old2new, new2old] = cleanMesh(mesh)
%CLEANMESH    Removes duplicate and unused points from mesh
%   [CLEANEDMESH, OLD2NEW, NEW2OLD] = ...
%                       CLEANMESH(MESH)
%   The arrays old2new and new2old are used for reindexing in the following
%   way:
%   cleanedmesh.elements == old2new(mesh.elements) 
%   mesh.elements == new2old(cleanedmesh.elements)
%       and
%   cleanedmesh.coordinates == mesh.coordinates(new2old,:)
%   Works for arbitrary-dimensional meshes.
%
% Author: Josef Kemetmueller - 16.12.13

% Compute remaining elements
remaining = unique(mesh.elements(:));
nRe = size(remaining,1);
nE = numElements(mesh);


old2new = zeros(nE,1);
new2old = zeros(nRe,1);

% Compute new elements
new2old(1:nRe) = remaining;
old2new(remaining) = 1:nRe;

mesh.elements = reshape(old2new(mesh.elements),nE,[]);
if isfield(mesh,'coordinates')
    mesh.coordinates = mesh.coordinates(new2old,:);
end
if isfield(mesh, 'elementgeneration')
    mesh = rmfield(mesh,'elementgeneration');
end
nBE = numBoundaryElements(mesh);
for j = 1:numBoundaries(mesh)
    mesh.bd(j).elements = reshape(old2new(mesh.bd(j).elements),nBE(j),[]);
end