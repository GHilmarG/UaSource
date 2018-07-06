function q = MeshQuality(coordinates,connectivity)

%  QUALITY: Approximate triangle quality. 
%
%  q = quality(p,t);
%
%  p: Nx2 array of nodal XY coordinates, [x1,y1; x2,y2; etc]
%  t: Mx3 array of triangles as indices, [n11,n12,n13; n21,n22,n23; etc]
%  q: Mx1 vector of triangle qualities. 0<=q<=1.

% Darren Engwirda - 2007.

% Nodes
[Nele,nod]=size(connectivity);

switch nod
    case 3
        p1 = coordinates(connectivity(:,1),:);
        p2 = coordinates(connectivity(:,2),:);
        p3 = coordinates(connectivity(:,3),:);
    case 6
        p1 = coordinates(connectivity(:,1),:);
        p2 = coordinates(connectivity(:,3),:);
        p3 = coordinates(connectivity(:,5),:);
    case 10
        p1 = coordinates(connectivity(:,1),:);
        p2 = coordinates(connectivity(:,4),:);
        p3 = coordinates(connectivity(:,7),:);
end

% Approximate quality
d12 = p2-p1;
d13 = p3-p1;
d23 = p3-p2;
q = 3.4641*abs(d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))./sum(d12.^2+d13.^2+d23.^2,2);

end      % quality()
