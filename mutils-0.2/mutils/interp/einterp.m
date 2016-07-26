function V = einterp(varargin)
%%EINTERP computes FEM interpolation at given points in an unstructured mesh.
%
% Usage: V_MARKERS = einterp(MESH, V, MARKERS, T)
%
% Input:
%  MESH        structure containing the mesh description
%   MESH.NODES mesh nodes coordinates (2 x nnod)
%   MESH.ELEMS elements definition (nnodel x nel)
%
%  V           nodal velocities (2 x nnod)
%  MARKERS     marker coordinates (2 x n_markers)
%  T           elements containing the markers
%
% Output:
%  V_MARKERES  interpolated values in MARKERS
%
% EINTERP uses SSE instructions, prefetching and OpenMP parallelization to
% achieve good performance on modern multi-core CPUs. Set the OMP_NUM_THREADS 
% variable to the desired number of CPU cores.
%
% Currently, the implementation is limited to 7-node triangular elements and
% two degrees of freedom per mesh node. More elements will be added as
% requested by the users.
%
% Elements supported:
%
%  7-node triangle, 2 dofs per node
%
%See also: TSEARCH, TSEARCH2

% Copyright 2012, Marcin Krotkiewski, University of Oslo

error ('MEX function not found');

end
