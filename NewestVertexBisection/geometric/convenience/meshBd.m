function mesh = meshBd(mesh, j)
%MESHBD    Returns a boundary as mesh structure
%   MESH = MESHBD(MESH, j) returns the j-th boundary and is basically a
%   short form for 
%   GENMESH(mesh.bd(j).elements, mesh.coordinates);
% 
%   See also:
%   GENMESH
%
%   Author: Josef Kemetmueller - 16.12.2013
if isfield(mesh, 'bd') && j<=length(mesh.bd)
    mesh = genMesh(mesh.bd(j).elements, mesh.coordinates);
else
    mesh = genMesh(zeros(0,dimMesh(mesh)), mesh.coordinates);
end