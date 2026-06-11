function [coordinates,connectivity]=RemoveDuplicateNodes(coordinates,connectivity,tolerance)
    
    % Removes duplicates nodes and duplicate elements and renumbers the connectivity 
    % so that it contains only numbers in the range=1:Nnodes 
    %
    % [coordinates,connectivity]=RemoveDuplicateNodes(coordinates,connectivity,tolerance)
    % 
    
    if size(connectivity,1) == 1 
        return
    end
    
    if nargin < 3
        tolerance=100*eps;
    end
    
    x=coordinates(:,1);  y=coordinates(:,2); 
    ix=round(x/tolerance); iy=round(y/tolerance);
    [~,n,m]=unique([ix,iy],'rows');
    
    x=coordinates(n,1);
    y=coordinates(n,2);
    coordinates=[x(:) y(:)];
    connectivity=m(connectivity);
         
    % now remove duplicates in connectivity
    [~,n,~]=unique(sort(connectivity,2),'rows');
    
    connectivity=connectivity(n,:);
    
end

