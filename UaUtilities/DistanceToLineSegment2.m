function Ind = DistanceToLineSegment2(p, A, B,nTol,aTol)

%%
%
% Determines if a point p is close to a line segment defined by the end points A
% and B.
%
%
%
% There are two input formats:
%
%   Ind = DistanceToLineSegment(p, A, B,tolerance)
%
% where A and B define the end points of a line segment, and (the easier to use)
% format:
%
%   Ind = DistanceToLineSegment(p, coo , [],tolerance)
%
% where coo is a vector defining the end points of individual connected line segments.
%
% A point is considered to be on a line segment if:
%
% # the distance normal to the line segment B-A (vector notation) to the point p is less than nTol, and
% # if the distance along the line segment is less than the length of the line
% segment using the tolerance aTol
%
% Usefull to determine if point p is on or close to the line segment B-A
% and to find all nodal points along the boundary defined by the pairwise joining
% of points in A and B.
%
% Inputs:
%
% p  : Nx2 matrix of with the (x,y) cooridnates of points
%
% A , B : Mx2 matrix defining start and end points of M line segments.
%
% coo   : Mx2 matrix defining the (x,y) coordinates of connected line segments.
%
% nTol and aTol are distances. Points
%
% Outputs:
%
% Ind : all points in p that are within tolerance to any of the line segments.
%
% The AlongDist and NormDist are min distances along and normal to the orientation
% of any of the line segments.  (The AlongDist is normalized.)
%
%
% Resonably fast if number of elements in p much greate than in A and B, i.e if many more points in
% the array p than there are number of line segments (points in A and B).
%
% Example:
%
%   p=[3 0 ; 0 0 ] ; a=[1 0] ; b=[5 0]; nTol=0.1;
%   [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, a, b,nTol)
%
% Example:
% if xx and yy are vectors defining (x,y) locations, I can find all boundary nodes
% that are along the line segments joining (x,y) in the following manner:
%
%   [Ind,AlongDist,NormDist] = DistanceToLineSegment([x(Boundary.Nodes) y(Boundary.Nodes)],[xx(:) yy(:)],[],tolerance);
%
% Boundary.Nodes(Ind) now gives me the nodal numbers of all nodes along the boundary defined by the line segments
%
% Example:
%
%   p=[0 0 ; 0.5 0  ; 1 0 ; 1.5 0 ] ; coo=[0 0 ; 1 0 ] ; tolerance=0.1 ;
%   Ind = DistanceToLineSegment2(p, coo, [],tolerance)
%
% gives Ind=[1 ; 2 ; 3]
%
%  Example:
%
%    N=1000 ; M=2 ; p=100*rand(N,2) ; x=0:100 ; y=10*sin(2*pi*x/50)+50 ; coo=[x(:) y(:)]  ; tolerance=10 ;
%    Ind = DistanceToLineSegment2(p, coo, [],tolerance) ; 
%    figure ; plot(p(:,1),p(:,2),'.k')  ; hold on ; plot(coo(:,1), coo(:,2),'r')
%    plot(p(Ind,1),p(Ind,2),'or')
%
%
% See also RangeSearch,  NeighborhoodMatrix   

if nargin<3 || isempty(B)
    B=[A(2:end,1) A(2:end,2)];
    A=[A(1:end-1,1) A(1:end-1,2)] ;
end

if nargin<4 ; nTol=eps ; end;

if nargin<5
    aTol=nTol;
end

M=size(A,1);
N=size(p,1);
Ind=false(N,1) ; 


parfor J=1:M
    
    a=A(J,:) ; b=B(J,:) ;
    
    d=b-a;
    l=norm(d);
    
    e=p;
    e(:,1)=p(:,1)-a(1);
    e(:,2)=p(:,2)-a(2);
    
    dalong=e*d'/l      ;               % distance, if point p within the `transverse' strip to b-a, then 0<AlongDist<1
    dnorm=abs(e*[d(2) -d(1)]'/l)   ;   % absolute distance normal to b-a
    
    temp=dalong>=-aTol & dalong<=l+aTol &  dnorm<=nTol;
    Ind=Ind | temp ;

end


end
