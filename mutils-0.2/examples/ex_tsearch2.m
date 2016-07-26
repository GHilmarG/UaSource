% Copyright 2012, Marcin Krotkiewski, University of Oslo

clear

rand('seed', 1);

%% Mesh - generate Delaunay triangulation of random points
n_points              = 5e4;
points                = rand(2, n_points);

% First, create a new Delaunay triangulation.
% This part is done much faster with triangle. However, MATLAB
% pointLocation routine used later in comparison requires 
% a Delaunay triangulation created by DelaunayTri.
DT   = DelaunayTri(points');
trep = TriRep(DT.Triangulation, points(1,:)', points(2,:)');


%% Markers - search for elements containing regularily placed markers
n_markers = 10e6;
[X, Y]    = meshgrid(linspace(0.01,0.99,ceil(sqrt(n_markers))));
markers = [X(:) Y(:)]';
n_markers = length(markers);

% -------------------------
% tsearch2
% -------------------------
ELEMS = int32(DT.Triangulation');
WS = [];
WS.NEIGHBORS = trep.neighbors()';
WS.NEIGHBORS(isnan(WS.NEIGHBORS)) = 0;
WS.NEIGHBORS = int32(WS.NEIGHBORS);
WS.xmin = 0;
WS.xmax = 1;
WS.ymin = 0;
WS.ymax = 1;

t = tic;
setenv('OMP_NUM_THREADS', '1');
T1 = tsearch2(points, ELEMS, markers, WS);
display(['tsearch2 (sequential): ' num2str(toc(t))]);

% Only the point location part is parallelized, not creation of the
% quad-tree structure. Hence, speedups are best when the number of mesh
% elements is much smaller than the number of located markers.
t = tic;
setenv('OMP_NUM_THREADS', '2');
T2 = tsearch2(points, ELEMS, markers, WS);
display(['tsearch2 (parallel): ' num2str(toc(t))]);

% compare results, must be identical
ndiff = numel(unique(T1-T2))-1;
if ndiff
    display(['WARNING: Maps obtained by sequential and parallel tsearch2 do not match. ' ...
        num2str(ndiff) ' different results.']);
end

%% Compare time and results to MATLAB tsearch
% -------------------------
% tsearch - Obsolete MATLAB
% -------------------------
try
    t = tic;
    T3 = tsearch(points(1,:), points(2,:), double(ELEMS)', markers(1,:), markers(2,:));
    display(['tsearch: ' num2str(toc(t))]);
    
    % compare results - might differ for nodes lying on edges
    ndiff = numel(unique(T1-int32(T3)))-1;
    if ndiff
        display(['WARNING: Maps obtained by tsearch and tsearch2 do not match. ' ...
            num2str(ndiff) ' different results.']);
    end
catch
    display('WARNING: tsearch does not work or has been removed');
end

% -------------------------
% pointLocation - CGAL, new MATLAB version
% -------------------------
t = tic;
T4 = pointLocation(DT, markers');
display(['pointLocation: ' num2str(toc(t))]);

% compare results - might differ for nodes lying on edges
ndiff = numel(unique(T1-int32(T4')))-1;
if ndiff
    display(['WARNING: Maps obtained by tsearch2 and pointLocation do not match. ' ...
        num2str(ndiff) ' different results.']);
end
