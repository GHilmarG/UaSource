function mesh = orientMesh(mesh, orientation)
%ORIENTMESH    Orients the mesh positively or negatively
%   MESH = ORIENTMESH(MESH,'inwards') or
%   MESH = ORIENTMESH(MESH,'outwards') will orient the elements of the mesh
%   and their boundaries in an inwards/outwards facing manner.
%   The resulting mesh CANNOT be used for bisection anymore, as the
%   orientations are important for the bisection scheme. It would have to
%   run through the genBisectionMesh3D process again.
%
%   Works for n-dimensional meshes in R^n for arbitrary n.
%
%   See also:
%	ORIENTMESHWITHOUTBOUNDARY
%
% Author: Josef Kemetmueller - 16.12.13
nB = numBoundaries(mesh);
assert((dimSpace(mesh) == dimMesh(mesh)), ...
                    'Dimensions unsuitable for orientation property.');
mesh = orientMeshWithoutBoundary(mesh, orientation);
if nB>0
    bd = getBoundary(mesh,orientation);
    for j = 1:nB
        [~,I,~] = intersect(sort(bd,2),sort(mesh.bd(j).elements,2),'rows');
        mesh.bd(j).elements = bd(I,:);
    end
end

% Remove elementgeneration field to mark the mesh as not being bisectable
% anymore
if isfield(mesh, 'elementgeneration')
    mesh = rmfield(mesh,'elementgeneration');
end