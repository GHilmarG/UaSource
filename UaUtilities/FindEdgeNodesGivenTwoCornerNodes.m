
function iE=FindEdgeNodesGivenTwoCornerNodes(MUA,FB)

% Finds nodes along edges given free boundary that only includes corner nodes. 
%
%    iE=FindEdgeNodesGivenTwoCornerNodes(MUA,FB)
%
% Note: 
% Calculate FB using matlab freeBoundary routine, for example as:
%
%    triGR=CreateFEmeshTriRep(MUA.connectivity,MUA.coordinates);
%    FB=freeBoundary(triGR);
%
%

switch MUA.nod
    
    case 3
        
        iE=[];
        
    case 6
        
        
        M = AdjacentNodes(MUA.connectivity);
        T=M(FB(:,1),:) & M(FB(:,2),:);
        [i,j,~]=find(T) ;
        [~,I]=sort(i);
        iE=j(I);
        
    case 10
        
        error('FindEdgeNodesGivenTwoCorner:CaseNotImplemented','Not yet implemented for 10 node elements.')
        
end


end
