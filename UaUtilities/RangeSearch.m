function [iA,iB,I,D2]=RangeSearch(cooA,cooB,r,D2)

%
% Finds all points in a set A within a given distance of points in set B, and
% all points in set B within a given distance from set A.
%
%    [iA,iB,I,D2]=RangeSearch(cooA,cooB,r)
%
% Inputs:
%  cooA   : n x 2 array of x and y locations of points in the plane.
%  cooB   : m x 2 array of x and y locations of points in the plane.
%     r   : a given distance (pos scalar)
%    D2   : Optional inputs. For repreated calls give
%           with same cooA and CooB but different r, give D2 from a previous
%           call as an input.
%
% 
% Outputs:
% iA    : logical array of all points in set A that are within the distance r
%         from any point in set B.
%
% iB    : logical array of all points in set B that are within the ditance r
%         from any point in set A.
%
% I     : Neigbourhood matrix. I(i,j) is true if point i in set A is within the
%         distance r from point j in set B.
%
% D2    : square-distance matrix.
%
%  Note:  This is a very simple implementation and limited to sets with not too
%  larger number of points. Usually memory rather than speed is the limitation.
%  A reasonably fast for sets with less than 10 000 points or so, but this is very much
%  memory dependent.
% 
% See also: DistanceToLineSegement2, NeighborhoodMatrix


if nargin<4 || isempty(D2)
    
    xa=single(cooA(:,1)) ;  ya=single(cooA(:,2)) ;
    xb=single(cooB(:,1)) ;  yb=single(cooB(:,2)) ;
    
    [Xa,Xb]=ndgrid(xa,xb) ;
    [Ya,Yb]=ndgrid(ya,yb);
    
    
    D2=(Xa-Xb).^2+(Ya-Yb).^2;
end

I=D2<r^2;
iA=any(I,2);

if nargout>1
    iB=any(I,1);
end

end