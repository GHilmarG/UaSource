function [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, A, B,tolerance)
    
    %%
    %
    % Determines if a point p is close to a line segment defined by the end points A
    % and B.
    %
    % Also calculates normal and along distances from a point to a line segment.
    %
    %   
    % There are two input formats:
    %
    %   [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, A, B,tolerance)
    %
    % where A and B define the end points of a line segment, and
    %
    %   [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, coo , [],tolerance)
    %
    % where coo is a vector defining the end points of individual connected line segments.
    %
    % A point is considered to be on a line segment if:
    %
    % # the distance normal to the line segment B-A (vector notation) to the point p is less than tolerance, and
    % # if the distance along the line segment is less than the length of the line segment using same tolerance
    %
    % Usefull to determine if point p is on or close to the line segment B-A
    % and to find all nodal points along the boundary defined by the pairwise joining
    % of points in A and B.
    %
    % Inputs:
    % 
    % p  : Nx2 matrix of points
    %
    % A , B : Mx2 matrix defining start and end points of M line segments.
    %
    % coo   : Mx2 matrix defining the (x,y) coordinates of connected line segments.
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
    %   p=[3 0 ; 0 0 ] ; a=[1 0] ; b=[5 0]; tolerance=0.1;
    %   [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, a, b,tolerance)
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
    %   Ind = DistanceToLineSegment(p, coo, [],tolerance)
    %
    % gives Ind=[1 ; 2 ; 3] 
    %
    %
    %
    
    if nargin<3 || isempty(B)
        B=[A(2:end,1) A(2:end,2)];
        A=[A(1:end-1,1) A(1:end-1,2)] ;
    end
    
    if nargin<4 ; tolerance=eps ; end;
    
    M=size(A,1);
    N=size(p,1);
    AlongDist=zeros(N,1)+Inf;
    NormDist=zeros(N,1)+Inf;
    
    
    Ind=[];
    for J=1:M
        a=A(J,:) ; b=B(J,:) ;
        
        d=b-a;
        l=norm(d);
        
        f(1)=d(2) ; f(2)=-d(1);
        e=p;
        e(:,1)=p(:,1)-a(1);
        e(:,2)=p(:,2)-a(2);
        
        dalong=e*d'/l^2      ;  % normalized distance, if point p within the `transverse' strip to b-a, then 0<AlongDist<1
        dnorm=abs(e*f'/l)   ;   % absolute distance normal to b-a
        
        Ind=[Ind;find(dalong>=-tolerance & dalong<=1+tolerance &  dnorm<=tolerance)];
        
        AlongDist=min(AlongDist,dalong);
        NormDist=min(NormDist,dnorm);
        
        
    end
   
    Ind=unique(Ind);
    
end
