% Copyright 2012, Marcin Krotkiewski, University of Oslo

clear;

%% Unstructured triangular mesh
% setup domain
points   = [0 0; 1 0; 1 1; 0 1]';
segments = [1 2; 2 3; 3 4; 4 1]';

% triangle options
opts.element_type     = 'tri7';
opts.gen_neighbors    = 1;
opts.triangulate_poly = 1;
opts.min_angle        = 30;
opts.max_tri_area     = 0.001;

% triangle input
tristr.points         = points;
tristr.segments       = int32(segments);

% generate mesh
MESH = mtriangle(opts, tristr);

% show mesh
ncorners = 3;
nel = length(MESH.ELEMS);
X = reshape(MESH.NODES(1,MESH.ELEMS(1:ncorners,:)), ncorners, nel);
Y = reshape(MESH.NODES(2,MESH.ELEMS(1:ncorners,:)), ncorners, nel);
C = zeros(size(X));
figure(1); clf;
h = patch(X, Y, C);
axis square
