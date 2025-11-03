
function [coordinates,connectivity]=ChangeElementType(coordinates,connectivity,To)

%%
%
% Changes element type from 3, 6, or 10 node triangles to 3, 6 or 10 node triangles
%
%   [coordinates,connectivity]=ChangeElementType(coordinates,connectivity,To)
%
% To : number of nodes of desired element type [3,6,10]

[Nele,nod]=size(connectivity);

if ~ismember(To,[3 6 10])
    fprintf('Changing to %-i nod element not possible\n',To)
    error('ChangeElementType')
end

From=nod;

if From==To ; return ; end

fprintf('Changing element type from %i-node to %i-node.\n',From,To)


% if I have 6 or 10 node elements, change all to 3 node elements
switch From
    case 6
        [coordinates,connectivity]=tri6to3(coordinates,connectivity);
        
    case 10
        [coordinates,connectivity]=tri10to3(coordinates,connectivity);
        %fprintf(' tri10to3 \n')
end

% now all elements are 3 nod elements


% change 3 node elements into the desired element type

switch To
    case 6
        %fprintf(' tri3to6 \n')
        [coordinates,connectivity]=tri3to6(coordinates,connectivity);
    case 10
        %fprintf(' tri3to10 \n')
        [coordinates,connectivity]=tri3to10(coordinates,connectivity);
end


end
