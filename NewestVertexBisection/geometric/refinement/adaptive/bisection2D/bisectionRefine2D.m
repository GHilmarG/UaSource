function mesh = bisectionRefine2D(mesh, markedElements)
%BISECTIONREFINE2D    Refines a mesh via bisection
%   MESH = BISECTIONREFINE2D(MESH, MARKEDELEMENTS) refines the elements
%   MARKEDELEMENTS via recursive bisection.
%
%   Works for 2-dimensional meshes in R^n.
%
%   Use the function BISECTIONREFINE instead of this function directly.
%
%   See also:
%   BISECTIONREFINE
%
% Author: Josef Kemetmueller - 20.03.13
assert(isfield(mesh, 'N0'), 'Mesh not generated using genBisectionMesh.');

nC = numCoordinates(mesh);
nE = numElements(mesh);
nB = numBoundaries(mesh);
bd2edges = cell(1,nB);
[edge2nodes, element2edges, bd2edges{1:nB}] = getEdges(mesh);
nEdges = size(edge2nodes,1);
if ischar(markedElements)
    assert(strcmpi(markedElements,'all'), 'Unknown option %s', markedElements);
    markedElements = 1:nE;
end
%% Find out which edges to refine
markedEdges = false(nEdges,1);
markedEdges(element2edges(markedElements,2)) = true; % (1,3) is the refinement edge.
additionalElements = 1; % Elements where one of the edges is marked, but not the refinement edge.
while any(additionalElements)
    el2isMarkedEdge = reshape(markedEdges(element2edges),[],3);
    additionalElements = ( ( el2isMarkedEdge(:,1) | el2isMarkedEdge(:,3) ) ...
                          & ~el2isMarkedEdge(:,2) );
    markedEdges(element2edges(additionalElements,2)) = true;
end
%% Generate new coordinates
MPs = @(EL, C) 0.5*(C(EL(:,1),:)+C(EL(:,2),:));
mesh.coordinates = [mesh.coordinates; MPs(edge2nodes(markedEdges,:), mesh.coordinates)];
newCoNums = nC+(1:nnz(markedEdges))';
%% Bisection rules
none     = ~el2isMarkedEdge(:,1) & ~el2isMarkedEdge(:,2) & ~el2isMarkedEdge(:,3);
bisec_2_ = ~el2isMarkedEdge(:,1) &  el2isMarkedEdge(:,2) & ~el2isMarkedEdge(:,3);
bisec12_ =  el2isMarkedEdge(:,1) &  el2isMarkedEdge(:,2) & ~el2isMarkedEdge(:,3);
bisec_23 = ~el2isMarkedEdge(:,1) &  el2isMarkedEdge(:,2) &  el2isMarkedEdge(:,3);
bisec123 =  el2isMarkedEdge(:,1) &  el2isMarkedEdge(:,2) &  el2isMarkedEdge(:,3);

bisec_2_rule = @(NODE, EDGE) [NODE(:,1), EDGE(:,2), NODE(:,2); ...
                              NODE(:,3), EDGE(:,2), NODE(:,2)];
bisec12_rule = @(NODE, EDGE) [NODE(:,1), EDGE(:,1), EDGE(:,2); ...
                              NODE(:,2), EDGE(:,1), EDGE(:,2); ...
                              NODE(:,3), EDGE(:,2), NODE(:,2)];
bisec_23rule = @(NODE, EDGE) [NODE(:,1), EDGE(:,2), NODE(:,2); ...
                              NODE(:,3), EDGE(:,3), EDGE(:,2); ...
                              NODE(:,2), EDGE(:,3), EDGE(:,2)];
bisec123rule = @(NODE, EDGE) [NODE(:,1), EDGE(:,1), EDGE(:,2); ...
                              NODE(:,2), EDGE(:,1), EDGE(:,2); ...
                              NODE(:,3), EDGE(:,3), EDGE(:,2); ...
                              NODE(:,2), EDGE(:,3), EDGE(:,2)];

%% We put the bisected elements right next to each other, for implicit binary tree usage in coarsening
idx = ones(nE,1);
idx(bisec_2_) = 2;
idx(bisec12_) = 3;
idx(bisec_23) = 3;
idx(bisec123) = 4;
idx = [1;1+cumsum(idx)];
%*** Generate new elements
edge2newNodeNum = zeros(nEdges, 1);
edge2newNodeNum(markedEdges) = newCoNums;
newNodes = reshape(edge2newNodeNum(element2edges),[],3);

where = @(bis,num)  bsxfun(@plus, reshape(idx(bis),[],1), 0:num);
newElements = zeros(idx(end)-1,3);
newElements(where(none,0),    :) = mesh.elements(none,:);
newElements(where(bisec_2_,1),:) = bisec_2_rule(mesh.elements(bisec_2_,:), newNodes(bisec_2_,:));
newElements(where(bisec12_,2),:) = bisec12_rule(mesh.elements(bisec12_,:), newNodes(bisec12_,:));
newElements(where(bisec_23,2),:) = bisec_23rule(mesh.elements(bisec_23,:), newNodes(bisec_23,:));
newElements(where(bisec123,3),:) = bisec123rule(mesh.elements(bisec123,:), newNodes(bisec123,:));

mesh.elements = newElements;
%% Refine boundaries
nBE = numBoundaryElements(mesh);
for j = 1:nB
    boundary = mesh.bd(j).elements;
    if ~isempty(boundary)
        newNodes = edge2newNodeNum(bd2edges{j});
        bisec = newNodes~=0;
        if any(bisec)
            idx = ones(nBE(j),1);
            idx(bisec) = 2;
            idx = [1;1+cumsum(idx)];
            newBoundary = zeros(idx(end)-1,2);
            newBoundary(idx(~bisec),:) = boundary(~bisec,:);
            newBoundary(bsxfun(@plus, idx(bisec), 0:1), :) = ...
                                    [boundary(bisec,1), newNodes(bisec); ...
                                     newNodes(bisec), boundary(bisec,2)];
            mesh.bd(j).elements = newBoundary;
        end

    end
end
