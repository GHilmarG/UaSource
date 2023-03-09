function In=IsInBox(Box,x,y)

% Find points inside of a box where
%
%  Box= [xleft xright ydown yup]; 
% 
% A much more general approach is provided by inpoly2
%
% 

In=x>Box(1) & x <Box(2)  & y> Box(3) & y < Box(4) ;



end