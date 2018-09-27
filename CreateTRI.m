function [TRI,DT]=CreateTRI(MUA,DT)

% [TRI,DT]=CreateTRI(MUA,DT)
% creates a delaunay triangulation from the nodal points, and then eliminates
% triangles outside of the mesh boundary. (Note; The elimination only works if
% there is just one connected boundary...this needs to be improved)
% DT is an optional input
%
% Note: Unless a Delaunay triangulation is really needed, just use 
% tri=TriFE(connectivity) 
% which is a much faster way of getting simple 3 node triangles.
%
%
x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ;

if nargin<2 || isempty(DT)
    DT=delaunayTriangulation(x,y);
end


ic=incenter(DT);

% note: This will fail if there are several polygons
[cnInt,on] = inpoly(ic,[x(MUA.Boundary.EdgeCornerNodes) y(MUA.Boundary.EdgeCornerNodes)]);

TRI=DT.ConnectivityList(cnInt,:);


% figure ; triplot(DT)
% hold on ; PlotMuaBoundary([],MUA,'r')
% plot(ic(:,1),ic(:,2),'*g')
% hold on ; plot(ic(cnInt,1),ic(cnInt,2),'+k')


%TRI=DT.ConnectivityList;

end