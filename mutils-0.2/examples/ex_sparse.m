% Copyright 2012, Marcin Krotkiewski, University of Oslo

clear;

%% Create unstructured triangular mesh
% unstructured mesh
opts.max_tri_area = 0.00002;
opts.element_type = 'tri3';
opts.gen_edges = 0;
tristr.points = [...
    -1 1 1 -1;...
    -1 -1 1 1];
tristr.segments = int32([...
    1 2 3 4;...
    2 3 4 1]);
tristr.segmentmarkers = int32([1 2 3 4]);
MESH = mtriangle(opts, tristr);


%% Create sparse matrix using sparse_create
display('sparse_create');
ndof   = 1;
nnodel = 3;

% create general 'symbolic' non-zero connectivity structure
opts.n_node_dof = ndof;
Aelem = 1;
t = tic;
Ag = sparse_create(MESH.ELEMS, Aelem, opts);   
display(['general symbolic sparse matrix (1 dof per node): ' num2str(toc(t))]);

% create symmetric lower-triangular connectivity structure
opts.symmetric = 1;
opts.n_node_dof = ndof;
Aelem = 1;
ELEMS=MESH.ELEMS;
t = tic;
As = sparse_create(MESH.ELEMS, Aelem, opts);   
display(['symmetric symbolic sparse matrix (1 dof per node): ' num2str(toc(t))]);

% assemble matrices for 3 dofs per node
ndof = 3;
nelemdof = nnodel*ndof;
nel = length(MESH.ELEMS);

% assemble dummy element matrices into general sparse matrix
% dummy symmetric element matrix, the same for all elements
Aelem = rand(nelemdof,nelemdof);
Aelem = Aelem+Aelem';
Aelem = repmat(Aelem(:), 1, nel);
opts.symmetric = 0;
opts.n_node_dof = ndof;
t = tic;
Ag = sparse_create(MESH.ELEMS, Aelem, opts);   
display(['assemble general sparse matrix (3 dof per node): ' num2str(toc(t))]);

% assemble dummy element matrices into general sparse matrix
% dummy symmetric element matrix, the same for all elements
indx_j = repmat(1:nnodel*ndof,nnodel*ndof,1); indx_i = indx_j';
indx_i = tril(indx_i);
indx_j = tril(indx_j);
indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
indx_j = indx_j(:); indx_j = indx_j(indx_j>0);
Aelem_s = Aelem((indx_j-1)*nelemdof + indx_i,:);
opts.symmetric = 1;
t = tic;
As = sparse_create(MESH.ELEMS, Aelem_s, opts);   
display(['assemble symmetric sparse matrix (3 dof per node): ' num2str(toc(t))]);

clear As Aelem_s;

%% standard MATLAB or sparse2 from SuiteSparse, if available
fprintf(1,'\n');
display('MATLAB version');

% prepare triplet indices
t = tic;

ELEMS  = MESH.ELEMS;
nnodel = size(ELEMS,1);
nel    = size(ELEMS,2);
nedof  = nnodel*ndof;
ELEM_DOF = zeros(nedof, nel,'int32');
for dof=1:ndof
    ELEM_DOF(dof:ndof:end,:) = ndof*(ELEMS-1)+dof;
end

% element dof indices
indx_j = repmat(1:nnodel*ndof,nnodel*ndof,1); indx_i = indx_j';
indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
indx_j = indx_j(:); indx_j = indx_j(indx_j>0);

A_i = ELEM_DOF(indx_i,:);
A_j = ELEM_DOF(indx_j,:);
A_i = double(A_i(:));
A_j = double(A_j(:));
display(['triplet indices: ' num2str(toc(t))]);

t=tic;
A = sparse(A_i, A_j, Aelem);
display(['assemble general sparse matrix (sparse): ' num2str(toc(t))]);

Ad = A-Ag;
Ad = Ad(:);
if max(abs(Ad(:))) > 1e-14
    display('WARNING: sparse matrices created by sparse and sparse_create differ');
end
clear Ag Ad;

if exist(['sparse2.' mexext]) == 3
    t=tic;
    A2 = sparse2(A_i, A_j, Aelem);
    display(['assemble general sparse matrix (sparse2): ' num2str(toc(t))]);
    
    % compare results
    Ad = A-A2;
    Ad = Ad(:);
    if max(abs(Ad(:))) > 1e-14
        display('WARNING: sparse matrices created by sparse and sparse2 differ');
    end
end
