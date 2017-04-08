function mesh = orientMeshWithoutBoundary(mesh, orientation)
%ORIENTMESHWITHOUTBOUNDARY    Orients the mesh positively or negatively
%   MESH = ORIENTMESHWITHOUTBOUNDARY(MESH,'inwards') or
%   MESH = ORIENTMESHWITHOUTBOUNDARY(MESH,'outwards') will orient the
%   elements of the mesh - but NOT its boundaries - in an inwards/outwards
%   facing manner. The resulting mesh CANNOT be used for bisection anymore,
%   as the orientations are important for the bisection scheme. It would
%   have to run through the genBisectionMesh3D process
%   again.
%
%   Works for n-dimensional meshes in R^n for arbitrary n.
%
%   See also:
%	ORIENTMESH
%
% Author: Josef Kemetmueller - 16.12.13

switch orientation
    case 'inwards'
        isCorrect = @(X) X<0;
    case 'outwards'
        isCorrect = @(X) X>0;
    otherwise
        error(['Unknown orientation option. ', ...
              'Only valid: ''inwards'' or ''outwards''']);
end
assert((dimSpace(mesh) == dimMesh(mesh)), ...
                    'Dimensions unsuitable for orientation property.');

signedVolumes = getSignedElementVolumes(mesh);
orientcorrect = isCorrect(signedVolumes);
% Change orientation
mesh.elements(~orientcorrect,[1,2]) = mesh.elements(~orientcorrect,[2,1]);
% Remove elementgeneration field to mark the mesh as not being bisectable
% anymore
if isfield(mesh, 'elementgeneration')
    mesh = rmfield(mesh,'elementgeneration');
end