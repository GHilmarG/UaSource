function hyperfaces = simplexHyperfaces(dim,orientation)
%SIMPLEXHYPERFACES    returns the hyperfaces of a simplex.
%   HYPERFACES = SIMPLEXHYPERFACES(DIM) returns the dim+1 hyperfaces of a
%   dim-dimensional simplex [1:dim+1]. The ordering of the nodes is
%   ascending and the i-th hyperface does not contain the node i.
%   This means that the orientations of the hyperfaces alternate between
%   inwards/outwards.
%   To get outwards facing orientation use:
%   HYPERFACES = SIMPLEXHYPERFACES(DIM,'outwards')
%   HYPERFACES = SIMPLEXHYPERFACES(DIM,'inwards')
%   This alternates every other face's orientation.
%
%   This function is meant to fix an ordering for the nodes of a simplex,
%   so that all functions depending on this ordering remain consistent.
%
% Author: Josef Kemetmueller - 20.03.14
tmp = (ones(dim+1,1)*(1:dim+1))';
hyperfaces = reshape(tmp(~eye(dim+1)),dim,[])'; % Remove diagonal
% The above sorting is not arbitrary:
% Files depending on the use of exactly this ordering include:
%   genBisectionMesh3D
% So don't change the ordering, this could destroy the bisection algorithm.
% (Why would you even?)
if exist('orientation','var')
    if strcmpi(orientation, 'outwards')
        hyperfaces(2:2:end,:) = hyperfaces(2:2:end,end:-1:1); % Orient outwards.
    elseif strcmpi(orientation, 'inwards')
        hyperfaces(1:2:end,:) = hyperfaces(1:2:end,end:-1:1); % Orient inwards.
    else
        error('Unknown option %s',orientation);
    end
end
end
