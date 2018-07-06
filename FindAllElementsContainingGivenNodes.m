
function E=FindAllElementsContainingGivenNodes(connectivity,NodeList)


% E=FindAllElementsContainingGivenNodes(connectivity,NodeList)
% Finds elements containing nodes in NodeList
%
% Returns a cell array
% Example :
%    E=FindAllElementsContainingGivenNodes(connectivity,[13 23])
%    returns E{13} listing the elements containing node 13
%    and E{23} listing elements containing node 23
%
%
for k=1:numel(NodeList)
    
    [I,~]=find(connectivity==NodeList(k));
    E(NodeList(k))={I};
    
end

end