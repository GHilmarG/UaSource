function mesh = bisectionCoarsen(mesh, markedElements, varargin)
%BISECTIONCOARSEN    Coarsens a bisection mesh
%   MESH = BISECTIONCOARSEN(MESH, markedElements) coarsens the elements
%   markedElements.
%   Default behaviour is to only coarsen nodes, whose entire node patches
%   have been marked for coarsening.
%   You can use 
%   MESH = BISECTIONCOARSEN(MESH, markedElements, 'force') to coarsen a
%   node patch, even though only one of its elements has been marked.
%   This however does only perform a single coarsening step, and especially
%   only if this is possible.
%   If you want to mark all elements you can use
%   MESH = BISECTIONCOARSEN(MESH, 'all');
%   Repeated application of above function will finally yield the original
%   mesh (For 2D this only holds for initially refined meshes. See:
%   GENBISECTIONMESH for more information)
%
%   Works for 1 and 2 dimensional meshes in R^d and 3-dimensional meshes in
%   R^3 only.
%
%   See also:
%   GENBISECTIONMESH
%   BISECTIONREFINE
%
% Author: Josef Kemetmueller - 20.03.13

assert(dimMesh(mesh)<=3, 'Sorry, there is no bisection for d=%d',dimMesh(mesh));
switch dimMesh(mesh)
    case 1
        mesh = bisectionCoarsen1D(mesh, markedElements); % No forcing in 1D
    case 2
        mesh = bisectionCoarsen2D(mesh, markedElements, varargin{:});
    case 3
        mesh = bisectionCoarsen3D(mesh, markedElements, varargin{:});
end