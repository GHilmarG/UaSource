function A = sparse_create(varargin)
%SPARSE_CREATE creates a Finite Element sparse matrix for a given mesh.
%
% Usage:
%
%  A = sparse_create(ELEMS, [Aelems=1], [opts], [dof_map=[]])
%
% Input:
%  ELEMS          definition of elements (nnodel x nel)
%  Aelems         element matrices.
%                 1 - builds a symbolic nnz structure (1 in every non-zero
%                 entry)
%                 a single element matrix - use this if element matrices
%                 are the same for all elements
%                 one element matrix per element - use this if element
%                 matrices are different for all elements
%                 See examples below.
%  dof_map        supply a node/dof map, e.g. node permutation
%
%  opts           structure containing the following fields:
%
%    n_row_entries  average number of entries per row in the sparse matrix.
%                   Pass -1 to use the default value for recognized elements.
%                   Default -1;
%    n_node_dof     number of degrees of freedom per node. 
%                   Default 1.
%    symmetric      create symmetric (1) or general (0) sparse matrix.
%                   Default 0.
%    gen_map        generate a map that maps entries in Aelems to
%                   matrix non-zero entries.
%                   Default 0.
%
% Output:
%  A              sparse matrix
%
%   Examples:
%
%    connectivity graph
%
%     A = sparse_create(ELEMS);
%
%    symmetric sparse, different element matrix for all elements:
%
%      opts = [];
%      opts.symmetric = 1;
%      Aelem = ones(size(ELEMS,1)*(size(ELEMS,1)+1)/2, size(ELEMS,2));
%      A = sparse_create(ELEMS, Aelem, opts);
%
%    symmetric sparse, common element matrix for all elements:
%
%      opts = [];
%      opts.symmetric = 1;
%      Aelem = ones(size(ELEMS,1)*(size(ELEMS,1)+1)/2,1);
%      A = sparse_create(ELEMS, Aelem, opts);
%
%    general sparse different element matrix for all elements:
%
%      Aelem = ones(size(ELEMS,1)^2, size(ELEMS,2));
%      A = sparse_create(ELEMS, Aelem);
%
%    general sparse, common element matrix for all elements:
%
%      Aelem = ones(size(ELEMS,1)^2,1);
%      A = sparse_create(ELEMS, Aelem);

% Copyright 2012, Marcin Krotkiewski, University of Oslo

error ('MEX function not found');

end
