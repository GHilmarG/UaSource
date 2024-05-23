function [isInside,isOnBounday]=InsideOutside(xy,boundary)

%%
% Determines if xy is inside or outside of boundary
%
% Just a simple wrapper around inpoly2 to take care of the possibility of
% several boundaries within the array 'boundary' seperated by NaNs
%
%   xy            : n x 2 array
%   Boundary      : m x 2 Array
%
% Example
% 
%
%   xy=[ 0.5 0.5 ; 0.5 1.5] ;
%   boundary=[ 0 0 ; 1 0 ; 1 1 ; 0 1 ; 0 0 ; nan nan ; 0 2 ; 1 2 ; 1 3 ; 0 3 ; 0 2] ;
%   [isInside,isOnBounday]=InsideOutside(xy,boundary) ; 
%
%   figure ; 
%   plot(boundary(:,1),boundary(:,2),LineWidth=2) ; hold on ; 
%   plot(xy(:,1),xy(:,2),"*r")
%   plot(xy(isInside,1),xy(isInside,2),"ok") ; axis padded
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