

function [coordinates,connectivity]=RemoveNodesNotInConnectivity(coordinates,connectivity)

%
% [coordinates,connectivity]=RemoveNodesNotInConnectivity(coordinates,connectivity)
% gets rid of nodes that are not in connectivity
% relabels all nodes so that numbers start with 1 with no gaps


if numel(coordinates)==0 
    
    warning('Ua:RemoveNodesNotInConnectivity:NoCoordinates','Number of coordinates is zero!')
    return 

end

if numel(connectivity)==0 
    
    warning('Ua:RemoveNodesNotInConnectivity:NoElements','Number of elements is zero!')
    return 

end


MeshNodes=unique(connectivity);  % MeshNodes are all the nodes in connectivity list

% create a sorted list of nodes that are not in connectivity
NodeLabels=1:length(coordinates);
ExtraNodes=NodeLabels; ExtraNodes(MeshNodes)=NaN;   ExtraNodes=sort(ExtraNodes); ExtraNodes(isnan(ExtraNodes))=[];

% I make my life easy here by using an existing m file to solve this and
% set the coordinates of all nodes that are not in connectivity equal to some node in connectivity
% and then get rid of duplicates

xtemp=coordinates(MeshNodes(1),1); ytemp=coordinates(MeshNodes(1),2);

coordinates(ExtraNodes,1)=xtemp; coordinates(ExtraNodes,2)=ytemp;
[coordinates,connectivity]=RemoveDuplicateNodes(coordinates,connectivity,1000*eps);

% PlotFEmesh(coordinates,connectivity)

end



