function mesh = genMesh(elements, coordinates, varargin)
%GENMESH    Constructor for a mesh structure
%   MESH = GENMESH(ELEMENTS, COORDINATES, BOUNDARY1, ...) constructs a
%   simplicial mesh, which is a structure containing the arrays ELEMENTS
%   and COORDINATES (and optionally a partioning of its boundary BD). These
%   represent the mesh by standard simplex-vertex format. The array
%   elements is given as numElements-by-dimMesh(mesh)+1 matrix, whereas
%   coordinates is given as numCoordinates-by-dimSpace matrix. Multiple
%   boundaries can be given. If no boundary is given, the call does not
%   generate a boundary. MESH = GENMESH(ELEMENTS, COORDINATES, 'any'), MESH
%   = GENMESH(ELEMENTS, COORDINATES, 'outwards') or MESH =
%   GENMESH(ELEMENTS, COORDINATES, 'inwards') generates
%   anywhere/inwards/outwards facing boundary facets via GETBOUNDARY. 
%
%   The following examples generate kuhn simplices, where we divide the
%   boundary into two parts, the first one being given explicitly, the
%   second one being computed implicitly. (This will result in a warning)
%       mesh = genMesh(1:4, [0,0,0;eye(3)], [1,2,3]);
%       OR
%       n = 3;
%       mesh = genMesh(1:n+1, [zeros(1,n);eye(n)], 1:n);       
%
%   Works for arbitrary-dimensional meshes. 
%   Orientation only works for n-dimensional meshes in R^n.
%
%   See also:
%   GENBISECTIONMESH
%   MESHBD
%   GETBOUNDARY
%
% Author: Josef Kemetmueller - 16.12.13
mesh = struct('elements', elements, ...
              'coordinates', coordinates);
if nargin==2 % When no boundary is given, there's nothing to do.
    return;
end


if ischar(varargin{1})
    if strcmpi(varargin{1},'any')
        varargin{1} = getBoundary(mesh);
    else
        varargin{1} = getBoundary(mesh,varargin{1});
    end
end
assert(all(cellfun(@(X) size(X,2)-1, varargin) == dimMesh(mesh)-1 | ...
           cellfun(@(X) size(X,2), varargin) == 0), ...
       'Boundaries must be of dimension dimMesh(mesh)-1.');
mesh.bd = cell2struct(varargin,'elements');

%% Now we check that everything is correct.
wholeBoundary = sort(getBoundary(mesh),2,'ascend');
givenBoundary = sort(cell2mat(reshape(varargin,[],1)),2,'ascend');
ngBE = size(givenBoundary,1);
% The next reshape is needed in case of empty givenBoundary
remainingBoundary = setdiff(wholeBoundary, ...
                            reshape(givenBoundary,[],dimMesh(mesh)),'rows');
for j = 1:numBoundaries(mesh)
   if isempty(mesh.bd(j).elements)
      mesh.bd(j).elements = zeros(0,dimMesh(mesh)); 
   end
end
assert(size(unique(givenBoundary,'rows'),1)==ngBE, ...
       'Some boundary facets are used multiple times.');
assert(isempty(setdiff(reshape(givenBoundary,[],dimMesh(mesh)),...
                       wholeBoundary,'rows')), ...
       'Given boundary is not subset of topological boundary.');
%% In case of remaining boundary facets we add them to the struct
if ~isempty(remainingBoundary)
    warning('Bis3D:NotWholeBoundary','Not all boundary facets were given.');
    mesh.bd(end+1).elements = remainingBoundary;
end

end