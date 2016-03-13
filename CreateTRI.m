function [TRI,DT]=CreateTRI(MUA,DT)

% [TRI,DT]=CreateTRI(MUA,DT)
% creates a delaunay triangulation from the nodal points, and then eliminates
% triangles outside of the mesh boundary. 
% DT is an optional input

x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ; 

if nargin<2 || isempty(DT)
    DT=delaunayTriangulation(x,y);
end

ic=incenter(DT);
[cnInt,on] = inpoly(ic,[x(MUA.Boundary.EdgeCornerNodes) y(MUA.Boundary.EdgeCornerNodes)]);
TRI=DT.ConnectivityList(cnInt,:);

end