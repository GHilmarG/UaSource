function [TRI,DT]=CreateTRI(MUA,DT)

% [TRI,DT]=CreateTRI(MUA,DT)
% creates a delaunay triangulation from the nodal points, and takes out
% those triangles within the convex hull but outside of the FE mesh
% Usefull, for example, when plotting with trisurf
% DT is optional

x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ; 

if nargin<2 || isempty(DT)
    DT=delaunayTriangulation(x,y);
end

ic=incenter(DT);
[cnInt,on] = inpoly(ic,[x(MUA.Boundary.EdgeCornerNodes) y(MUA.Boundary.EdgeCornerNodes)]);
TRI=DT.ConnectivityList(cnInt,:);

end