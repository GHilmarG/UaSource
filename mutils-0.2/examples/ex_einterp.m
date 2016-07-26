function ex_einterp

% Copyright 2012, Marcin Krotkiewski, University of Oslo

%% Generate unstructured mesh
opts.max_tri_area = 0.0001;
opts.element_type = 'tri7';
opts.gen_edges = 0;
opts.gen_neighbors = 1;
tristr.points = [...
    -2 2 2 -2;...
    -1 -1 1 1];
tristr.segments = int32([...
    1 2 3 4;...
    2 3 4 1]);

MESH = mtriangle(opts, tristr);

if size(MESH.ELEMS,1)==6
    % triangle does not create 7-node elements.
    % Add central nodes to the element structure, number last
    nel = length(MESH.ELEMS);
    ncorners = 3;
    MESH.ELEMS(end+1,:)  = max(MESH.ELEMS(:))+int32([1:nel]);
    MESH.NODES = [MESH.NODES [...
        mean(reshape(MESH.NODES(1, MESH.ELEMS(1:ncorners,:)), ncorners, nel));...
        mean(reshape(MESH.NODES(2, MESH.ELEMS(1:ncorners,:)), ncorners, nel))]];
end

%% Markers - generate markers, search for elements containing the markers
n_markers = 5e6;

% lexigraphic marker numbering
[X, Y]    = meshgrid(linspace(0,1,ceil(sqrt(n_markers))));
markers = [X(:) Y(:)]';

% scale the markers coordinates to fit the domain
markers(1,:)   = 4*markers(1,:)-2;
markers(2,:)   = 2*markers(2,:)-1;
n_markers      = length(markers);

% locate markers in elements using tsearch2
WS = [];
WS.NEIGHBORS = MESH.NEIGHBORS;
WS.xmin = -2;
WS.xmax =  2;
WS.ymin = -1;
WS.ymax =  1;

setenv('OMP_NUM_THREADS', '1');
t=tic;
T = tsearch2(MESH.NODES, MESH.ELEMS(1:3, :), markers, WS);
display(['tsearch2: ', num2str(toc(t))]);


%% FEM interpolation of 2D velocities in markers
% generate random velocity field
V = 1+rand(size(MESH.NODES));

% einterp MEX function
t=tic;
Vm_seq = einterp(MESH, V, markers, T);
display(['einterp MEX (sequential): ', num2str(toc(t))]);

setenv('OMP_NUM_THREADS', '2');
t=tic;
Vm = einterp(MESH, V, markers, T);
display(['einterp MEX (parallel): ', num2str(toc(t))]);

if unique(Vm-Vm_seq) ~= 0
    %save test_einterp;
    merror('sequential and parallel einterp results differ');
end

% native MATLAB
t=tic;
eta  = local_coordinates(MESH,markers,T);
eta1 = eta(1,:);
eta2 = eta(2,:);
eta3 = 1-eta1-eta2;

eta1eta2eta3 = eta1.*eta2.*eta3;
N = zeros(n_markers,7);
N(:,1) = eta1.*(2*eta1-1) + 3*eta1eta2eta3;
N(:,2) = eta2.*(2*eta2-1) + 3*eta1eta2eta3;
N(:,3) = eta3.*(2*eta3-1) + 3*eta1eta2eta3;
N(:,4) = 4*eta2.*eta3 - 12*eta1eta2eta3;
N(:,5) = 4*eta1.*eta3 - 12*eta1eta2eta3;
N(:,6) = 4*eta1.*eta2 - 12*eta1eta2eta3;
N(:,7) =                27*eta1eta2eta3;

ELEMS = MESH.ELEMS(:,T);
Vx = V(1,:);
Vy = V(2,:);
vx = sum(Vx(ELEMS).*N');
vy = sum(Vy(ELEMS).*N');
display(['einterp MATLAB: ', num2str(toc(t))]);


%% compare results: MEX vs MATLAB
Vs = [vx; vy];
display(['Maximum difference between MATLAB and MEX implementations: ' num2str(max(abs(unique(Vm-Vs))))]);


%% functions
    % function that computes local element coordinates of randomly placed
    % markers in an unstructured triangular mesh.
    function eta = local_coordinates(MESH,points,point_elems)
        ndim = size(MESH.NODES, 1);
        nnod = size(MESH.NODES, 2);
        nel  = length(MESH.ELEMS);
        
        ENOD_X = reshape(MESH.NODES(1,MESH.ELEMS(1:3,:)), 3,nel);
        ENOD_Y = reshape(MESH.NODES(2,MESH.ELEMS(1:3,:)), 3,nel);
        
        area  = ENOD_X(2,:).*ENOD_Y(3,:) - ENOD_X(3,:).*ENOD_Y(2,:) + ...
            ENOD_X(3,:).*ENOD_Y(1,:) - ENOD_X(1,:).*ENOD_Y(3,:) + ...
            ENOD_X(1,:).*ENOD_Y(2,:) - ENOD_X(2,:).*ENOD_Y(1,:);
        
        ENOD_X_LONG = ENOD_X(:,point_elems);
        ENOD_Y_LONG = ENOD_Y(:,point_elems);
        
        eta = zeros(ndim,length(points));
        eta(1,:)  = ENOD_X_LONG(2,:).*ENOD_Y_LONG(3,:) - ENOD_X_LONG(3,:).*ENOD_Y_LONG(2,:) + ...
            ENOD_X_LONG(3,:).*points(2,:) - points(1,:).*ENOD_Y_LONG(3,:) + ...
            points(1,:).*ENOD_Y_LONG(2,:) - ENOD_X_LONG(2,:).*points(2,:);
        
        eta(2,:)  = ENOD_X_LONG(3,:).*ENOD_Y_LONG(1,:) - ENOD_X_LONG(1,:).*ENOD_Y_LONG(3,:) + ...
            ENOD_X_LONG(1,:).*points(2,:) - points(1,:).*ENOD_Y_LONG(1,:) + ...
            points(1,:).*ENOD_Y_LONG(3,:) - ENOD_X_LONG(3,:).*points(2,:);
        
        area_long = area(point_elems);
        
        eta(1,:) = eta(1,:)./area_long;
        eta(2,:) = eta(2,:)./area_long;
    end
end