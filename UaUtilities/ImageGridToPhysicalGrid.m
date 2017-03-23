function [x,y,Z]=ImageGridToPhysicalGrid(x,y,Z)

% Often an image Z it is alinged such that 
%
% Z(i,j)  is  at (x(j),y(i)) with y values decreasing and x values increasing
% 
% This function flips to 
%
% Z(i,j) at (x(i),y(j)) with x and y increasing.
%
% Contouring will require transposing Z, i.e. contour(x,y,Z') 
%
%

x=x(:); 
y=y(:);

y=flipud(y);
Z=flipud(Z);

Z=Z';


end