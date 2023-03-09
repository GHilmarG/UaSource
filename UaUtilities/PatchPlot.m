function [PatchObject,ColorbarHandel,tri]=PatchPlot(x,y,z,varargin)

%% 
DT=delaunayTriangulation(x,y);
tri=DT.ConnectivityList ;


PatchObject=patch('faces',tri,'vertices',[x y],...
    'FaceVertexCData',z,'CDataMapping','scaled','EdgeColor','none','FaceColor','interp',varargin{:}) ;

axis equal
ColorbarHandel=colorbar;



return
