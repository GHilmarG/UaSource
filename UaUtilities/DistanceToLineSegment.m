function [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, A, B,tolerance)
    
    
    % determines if both
    % 1) the distance normal to the line segment B-A (vector notation) to the point p is less than tolerance, and
    % 2) if the distance along the line segment is less than the length of the line segment using same tolerance
    %
    % [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, A, B,tolerance)
    % [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, coo , [],tolerance)
    %
    % usefull to determine if point p is on or close to the line segment B-A
    % and to find all nodal points along the boundary defined by the pairwise joining
    % of points in A and B.
    %
    %
    % A and B are arrays of points in the plane, p a nx2 matrix of xy locations
    % returns I as a index array
    %
    % p  : Nx2 matrix of points
    % A , B ; Mx2 matrix defining start and end points of M line segments
    %
    % [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, coo , [],tolerance
    %  finds p along coo with tolerance
    %
    %
    % Ind : all points in p that are close to some of the line segments defined by A and B
    %
    %  Resonably fast if number of elements in p much greate than in A and B, i.e if many more points in 
    %  the array p than there are number of line segments (points in A and B).
    %
    % example:
    %         p=[3 0 ; 0 0 ] ; a=[1 0] ; b=[5 0]; tolerance=0.1;
    %         [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, a, b,tolerance)
    %
    % example: 
    % if xx and yy are vectors defining (x,y) locations, I can find all boundary nodes 
    % that are along the line segments joining (x,y) in the following manner:
    %
    % A=[xx(1:end-1) yy(1:end-1)] ; B=[xx(2:end) yy(2:end)]; tolerance=eps;
    % [Ind,AlongDist,NormDist] = DistanceToLineSegment([x(Boundary.Nodes) y(Boundary.Nodes)],[xx(:) yy(:)],[],tolerance);
    % Boundary.Nodes(Ind) now gives me the nodal numbers of all nodes along the boundary defined by the line segments 
    % 
    %
    
    
    if nargin<3 || isempty(B)
        B=[A(2:end,1) A(2:end,2)];
        A=[A(1:end-1,1) A(1:end-1,2)] ;
    end
    
    if nargin<4 ; tolerance=eps ; end;
    
    N=size(A,1);

    Ind=[];
    for J=1:N
        a=A(J,:) ; b=B(J,:) ;
        
        d=b-a;
        l=norm(d);
        
        f(1)=d(2) ; f(2)=-d(1);
        e=p;
        e(:,1)=p(:,1)-a(1);
        e(:,2)=p(:,2)-a(2);
        
        AlongDist=e*d'/l^2      ;  % normalized distance, if point p within the `transverse' strip to b-a, then 0<AlongDist<1
        NormDist=abs(e*f'/l)   ;   % absolute distance normal to b-a
        % is 0<= projection <= l ?
        
        Ind=[Ind;find(AlongDist>=-tolerance & AlongDist<=1+tolerance &  NormDist<=tolerance)];
        
    end
   
    Ind=unique(Ind);
    
end
