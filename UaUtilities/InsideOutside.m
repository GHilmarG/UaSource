function [isInside,isOnBounday]=InsideOutside(xy,boundary)

%%
% Determines if xy is inside or outside of boundary
%
% Just a simple wrapper around inpoly2 to take care of the possibility of
% several boundaries within the array 'boundary' seperated by NaNs
%
%%

Kisnan=find(isnan(boundary(:,1)));


if numel(Kisnan)==0
    
    [isInside,isOnBounday]=inpoly2(xy,boundary);
    
else
    
    isInside=false(size(xy(:,1)));
    isOnBounday=false(size(xy(:,1)));
    i1=1;
    
    for j=1:numel(Kisnan)
        
        i2=Kisnan(j)-1;
        [J,K]=inpoly2(xy,boundary(i1:i2,:));
        i1=i2+2;
        
        isInside=isInside | J;
        isOnBounday=isOnBounday | K;
    end
    
    [J,K]=inpoly2(xy,boundary(i1:end,:));
    isInside=isInside | J;
    isOnBounday=isOnBounday | K;
end


end