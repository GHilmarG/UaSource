function mesh = bisectionRefine1D(mesh, markedElements)
%BISECTIONREFINE1D    Refines a mesh via bisection
%   MESH = BISECTIONREFINE1D(MESH, MARKEDELEMENTS) refines the elements
%   MARKEDELEMENTS via recursive bisection.
%
%   Works for 1-dimensional meshes in R^n.
%
%   Use the function BISECTIONREFINE instead of this function directly.
%
%   See also:
%   BISECTIONREFINE
%
% Author: Josef Kemetmueller - 20.03.13

%% Well, THAT's easy!
warning('1D Refinement is programmed in a way, that neighboring elements can have huge differences in size. TODO...');

assert(isfield(mesh, 'N0'), 'Mesh not generated using genBisectionMesh.');
nC = numCoordinates(mesh);
nE = numElements(mesh);
if ischar(markedElements)
    assert(strcmpi(markedElements,'all'), 'Unknown option %s', markedElements);
    markedElements = 1:nE;
end
if islogical(markedElements)
    assert(length(markedElements)==numElements(mesh),'markedElements not in correct range');
    nMarked = nnz(markedElements);
else
    assert(~isempty(markedElements) && all(markedElements>=0) && all(markedElements<=numElements(mesh)),...
       'markedElements must be in range of number of elements.');
    nMarked = length(markedElements);
end

MP = @(EL,C) 0.5*(C(EL(:,1),:)+C(EL(:,2),:));
newElNum = nE+(1:nMarked)';
newCoNum = nC+(1:nMarked)';
mesh.coordinates(newCoNum,:) = MP(mesh.elements(markedElements,:), ...
                                  mesh.coordinates);
mesh.elements([markedElements(:);newElNum],:) = ...
                        [mesh.elements(markedElements,1), newCoNum;
                         newCoNum, mesh.elements(markedElements,2)];                     
%% 1D boundaries don't change! YIPPIEH!