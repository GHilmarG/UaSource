% Copyright 2012, Marcin Krotkiewski, University of Oslo

clear;

%% Unstructured triangular mesh
% setup domain
points   = [0 0; 1 0; 1 1; 0 1]';
segments = [1 2; 2 3; 3 4; 4 1]';

% add a circular inclusion
no_pts_incl  = 150;
radius       = 0.1;
x_min        = 0; 
x_max        = 1;
y_min        = 0; 
y_max        = 1;
alpha        = 0.5;

theta        = linspace(0,2*pi,no_pts_incl);
theta(end)   = [];
xx           = cos(theta);
yy           = sin(theta);
center_x     = alpha*x_max+(1-alpha)*x_min;
center_y     = 0.5*(y_max+y_min);
INCLUSION    = [center_x + radius*xx; center_y + radius*yy];
no_pts       = size(INCLUSION,2);
pts_u        = 4 + no_pts;
pts_l        = 5;
INCLUSION_s  = [pts_l:pts_u;pts_l+1:pts_u+1];
INCLUSION_s(2,end) = pts_l;
points       = [points INCLUSION];
segments     = [segments INCLUSION_s];

% triangle options
opts.element_type     = 'tri3';
opts.gen_neighbors    = 1;
opts.triangulate_poly = 1;
opts.min_angle        = 30;
opts.max_tri_area     = 0.01;

% triangle input
tristr.points         = points;
tristr.segments       = int32(segments);

% generate mesh
MESH = mtriangle(opts, tristr);

% show mesh
ncorners = 3;
nel = length(MESH.ELEMS);
X   = reshape(MESH.NODES(1,MESH.ELEMS(1:ncorners,:)), ncorners, nel);
Y   = reshape(MESH.NODES(2,MESH.ELEMS(1:ncorners,:)), ncorners, nel);
C   = zeros(size(X));
figure(1); clf;
h = patch(X, Y, C);
axis square


%% Build quad-tree based on mesh nodes
qtree = quadtree('create', MESH.NODES, x_min, x_max, y_min, y_max, 5);

% save qtree in a VTK file
quadtree('vtkwrite', qtree, 'example');


%% Renumber random cloud of points using Z-curve ordering
% randomize markers
n_markers = 1e6;
markers   = 0.01+0.98*rand(2, n_markers);

% Create a quadtree based on the markers.
% The quadtree is only roughly adapted (large number of markers
% in a quadrant allowed). The quadtree is used to reorder the
% markers to achieve 2D spatial locality.
qtree_markers = quadtree('create', markers, x_min, x_max, y_min, y_max, 4096);
    
% find marker ordering based on quadtree traversal (Z-curve)
I = quadtree('reorder', qtree_markers);
markers = markers(:,I);
