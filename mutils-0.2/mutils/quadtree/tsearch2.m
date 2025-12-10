function [T, WS, stats] = tsearch2(NODES, TRI, POINTS, WS, T)
%%TSEARCH2 locates points in triangles using quadtree. 
%
% Usage:
% 
%  [T, WS, stats] = tsearch2(NODES, TRI, POINTS, WS, [T])
%
% Input: as in tsearch, except:
%  
%   TRI            triangle definitions (type int32). 
%                  TRI needs not be a Delaunay triangulation. 
%
%   WS is a structure containing the following fields:
%
%    WS.NEIGHBORS   a 3 X nel array of triangle neigbors. Order of
%                   neighbors matters: neighbor 1 is opposite node 1, etc.
%                   If there is no corresponding neighbor, NEIGHBORS should 
%                   contain 0. See HELP QUADTREE for details. Obligatory.
%    WS.ELEMS_CENTERS coordinates of the central points of the triangles. Optional.
%    WS.xmin
%    WS.xmax
%    WS.ymin
%    WS.ymax        Spatial domain extents of the MESH. Optional. 
%                   Default from 0 to 1.
%    WS.qtree       Quad-tree structure created by an earlier call to
%                   TSEARCH2. qtree can be reused between searches done on 
%                   the same mesh, which may speed things up depending on the
%                   problem size. Optional.
%
%   T               Approximate/incomplete point locations (type int32). 
%                   Should contain 0 if approximate point location is 
%                   not known. Optional.
% 
% Output:
%
%   T               Triangle IDs. 0 if no triangle contains a given point.
%   WS              Updated workspace, as in input.
%   stats           Useful point location statistics, e.g. number of
%                   point-in-triangle tests.
%
%   Note taht the input/output types are different than in tsearch: 
%   T and TRI are int32.
%
% TSEARCH2 uses quad-tree structure to quickly locate the triangle
% containing a given point. First, qtree structure is built for the centers
% of triangular elements. Next, quadtree('locate', ...) is invoked.
%
% Parallelized on SMP machines using OpenMP. Set the OMP_NUM_THREADS 
% variable to the desired number of CPU cores.
%
% See also: QUADTREE, TSEARCH

% Copyright 2012, Marcin Krotkiewski, University of Oslo

%% Check number of parameters, their types and sizes
% Minimum and maximum number of parameters
error(nargchk(4, 5, nargin, 'struct'))

% Optional parameters - check if supplied, set to [] if not.
if nargin < 5;  T  = int32([]); end

% Check types of all parameters. Syntax similar to validateattributes
validateattributes(NODES,  {'double'}, {'size', [2 NaN]});
validateattributes(TRI,    {'int32'},  {'size', [3 NaN]});
validateattributes(POINTS, {'double'}, {'size', [2 NaN]});

if ~isfield(WS, 'NEIGHBORS')
    error('WS must contain NEIGHBORS field.');
end
if ~isfield(WS, 'xmin');  WS.xmin = 0; end
if ~isfield(WS, 'xmax');  WS.xmax = 1; end
if ~isfield(WS, 'ymin');  WS.ymin = 0; end
if ~isfield(WS, 'ymax');  WS.ymax = 1; end
if ~isfield(WS, 'qtree'); WS.qtree = []; end

if ~isempty(T)
    validateattributes(T, {'int32'}, {'vector'});
end


%% Work
% the below can be removed and is only here for tsearch compatibility
MESH.NODES = NODES;
MESH.ELEMS = TRI;
MESH.NEIGHBORS = WS.NEIGHBORS;

% Is a quadtree structure supplied?
nel = size(TRI, 2);
if ~isempty(WS.qtree);
    % verify that qtree was created for a system 
    % with the same number points
    if nel~=WS.qtree.n_points
        error('qtree data structure is inconsistent with passed triangulation.');
    end
else
    if ~isfield(WS, 'ELEMS_CENTERS')
        WS.ELEMS_CENTERS = [...
            mean(reshape(MESH.NODES(1, MESH.ELEMS), 3, nel));...
            mean(reshape(MESH.NODES(2, MESH.ELEMS), 3, nel))];
    end
    WS.qtree = quadtree('create', WS.ELEMS_CENTERS, WS.xmin, WS.xmax, WS.ymin, WS.ymax, 2);
end

% location with the help of quadtree and element neigobors
[T, stats] = quadtree('locate', WS.qtree, MESH, POINTS, T);

end
