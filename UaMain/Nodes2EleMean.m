function [MeanEleValues]=Nodes2EleMean(connectivity,NodalValues)
    
    %
    % [MeanEleValues]=Nodes2EleMean(connectivity,NodalValues)
    % Takes nodal values and returns average of nodal values for each element
    %
    
    [Nele,nod]=size(connectivity);
    MeanEleValues=mean(reshape(NodalValues(connectivity),Nele,nod),2);
    
end
