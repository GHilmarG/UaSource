function [coordinates,connectivity]=RemoveNodesNotPartOfAnyElement(coordinates,connectivity)

% Removes nodes that are not part of any element,  renumbers the nodes and modifies connectivity accordingly
% Note: if there are duplicate nodes in coordinates and these nodes then appear separatly in connectivity,
%       then these two nodes will still be in the modified coordinates and connectivity

% CtrlVar.PlotEleLabels=1;
% figure
% PlotFEmesh(coordinates,connectivity,CtrlVar);
% title('original')


b=unique(connectivity(:));
ind(b(1:numel(b)))=1:numel(b);

coordinates=coordinates(b,:);
connectivity=ind(connectivity);

% figure
% PlotFEmesh(coordinates,connectivity,CtrlVar);
% title('final')

end

