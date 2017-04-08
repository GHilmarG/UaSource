function normals = getNormals(mesh, orientation)
%GETNORMALS    Orthonormal vectors of each elements hyperfaces.
%   NORMALS = GETNORMALS(MESH, 'outwards') or
%   NORMALS = GETNORMALS(MESH, 'inwards') returns the outwards/inwards
%   facing normals of all elements hyperfaces. NORMALS is a (DIM+1)-cell
%   array, so that NORMALS{j}(i,:) is the normal of the face opposite to
%   the node MESH.ELEMENTS(i,j) on the element MESH.ELEMENTS(i,:).
%
%   Works for arbitrary-dimensional meshes.
%
% Author: Josef Kemetmueller - 20.03.14

hatGrads = getHatGrads(mesh);
lengths = @(X) sqrt(sum(X.^2,2));
normalize = @(X) bsxfun(@rdivide, X, lengths(X));

normals = cellfun(normalize, hatGrads, 'UniformOutput', false);
if strcmpi(orientation, 'inwards')
    return;
elseif strcmpi(orientation, 'outwards')
    normals = cellfun(@uminus, normals, 'UniformOutput', false);
else
    error('Unknown option %s',orientation);
end
