function [varargout] = getElementsOfBoundary(mesh)
%GETELEMENTSOFBOUNDARY   To which element are the boundary facets attached to
%   [BDONE2ELEMENT, BDTWO2ELEMENT, ...] = getElementsOfBoundary(mesh)
%   returns arrays BDX2ELEMENT containing indices describing to which
%   element a certain boundary facet belongs.
%   mesh.elements(BDONE2ELEMENT(i),:) is the element which is attached to
%   the boundary facet mesh.bd(1).elements(i,:);
%
%   See also:
%   GETHYPERFACES
%
% Author: Josef Kemetmueller - 01.04.13

bd2faces = cell(1, numBoundaries(mesh));
[face2nodes, element2faces, bd2faces{:}] = getHyperfaces(mesh);
face2element = zeros(size(face2nodes,1), 1);
face2element(element2faces) = repmat((1:numElements(mesh))',1,4);

varargout = cell(1, numBoundaries(mesh));
for j = 1:numBoundaries(mesh);
    varargout{j} = face2element(bd2faces{j});
end