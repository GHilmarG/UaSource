function mesh = bisectionRefine(mesh, markedElements)
%BISECTIONREFINE    Refines a mesh via bisection
%   MESH = BISECTIONREFINE(MESH, MARKEDELEMENTS) refines the elements
%   MARKEDELEMENTS via recursive bisection.
%   This requires a compatible mesh, see GENBISECTIONMESH.
%
%   Notes for the 3D case:
%   The orientations of the refined elements and boundary faces will NOT be
%   uniform. This is a limitation of the underlying scheme. The
%   ordering/orientation of the elements and the nodes MUST NOT be changed,
%   as this will mess up the refinement routine and will probably lead to
%   MATLAB crashing. If you need a positively oriented mesh, you can get an
%   oriented copy of the mesh using ORIENTMESH.
%
%   Works for 1 and 2 dimensional meshes in R^d and 3-dimensional meshes in
%   R^3 only.
%
%   See also:
%   GENBISECTIONMESH
%   BISECTIONREFINE
%
% Author: Josef Kemetmueller - 20.03.13

assert(dimMesh(mesh)<=3, 'Sorry, there is no bisection implemented for d=%d',dimMesh(mesh));
switch dimMesh(mesh)
    case 1
        mesh = bisectionRefine1D(mesh, markedElements);
    case 2
        mesh = bisectionRefine2D(mesh, markedElements);
    case 3
        mesh = bisectionRefine3D(mesh, markedElements);
end