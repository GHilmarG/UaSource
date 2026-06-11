function boundary = getBoundary(mesh, varargin)
%GETBOUNDARY    Returns the topological boundary of the mesh.
%   BOUNDARY = GETBOUNDARY(MESH) returns the topological boundary.
%
%   BOUNDARY = GETBOUNDARY(MESH, 'outwards') or
%   BOUNDARY = GETBOUNDARY(MESH, 'inwards') returns the boundary in a way
%   that the normals are oriented outwards resp. inwards.
%
%   Works for arbitrary-dimensional meshes. 
%   Orientation only works for n-dimensional meshes in R^n.
%
% Author: Josef Kemetmueller - 16.12.13
% Note: freeBoundary(TR) should do the same - but without orientation.

nE = numElements(mesh);

if (nE==0)
    boundary = [];
    return;
end

%% Correct orientation
if (nargin>1)
    mesh = orientMeshWithoutBoundary(mesh, varargin{1});
end

%% Compute hyperfaces of the mesh
faceorder = simplexHyperfaces(dimMesh(mesh),'outwards');
hyperfaces = reshape(mesh.elements(:,faceorder),[],dimMesh(mesh));
sortedfaces = sort(hyperfaces,2);

[~,I,J] = unique(sortedfaces,'rows');
boundary = hyperfaces(I(accumarray(J,1)==1),:);
