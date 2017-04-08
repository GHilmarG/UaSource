function patchVols = getPatchVolumes(mesh)
%GETPATCHVOLUMES    Volumes of the node patches.
%   PATCHVOLS = GETPATCHVOLUMES(MESH) returns the vector of the sums of the
%   volumes of the simplices adjacent to each node.
%
%   Works for arbitrary-dimensional meshes.
%
%   Author: Josef Kemetmueller - 16.12.2013
nE = numElements(mesh);
% Patch integrals of the constant 1-function.
patchVols = integrateP0(mesh,ones(nE,1),'patchwise');