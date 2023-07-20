function [coordinates,connectivity,p]=NodalSweep(coordinates,connectivity,alpha)
        
% p is the permutation of nodal numbers 
%
% p provides the mapping between old and new nodal numbers
%
% old node number i is now node p(i)

    temp=coordinates(:,1)*cos(alpha)+coordinates(:,2)*sin(alpha);
    [~,p]=sort(temp);
    
    p2=p*0; p2(p)=1:length(p2);
    coordinates(:,2)=coordinates(p,2);
    coordinates(:,1)=coordinates(p,1);
    connectivity=p2(connectivity);
            
    
end

