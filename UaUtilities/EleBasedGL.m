
function [xGLele,yGLele,triGR,FB,xEdge,yEdge,nx,ny,ds]=EleBasedGL(CtrlVar,MUA,GF,DoPlots)

%%
% Finds a grounding line defined as the free boundary of all grounded elements.
%
%    [xGLele,yGLele,triGR,FB,xEdge,yEdge,nx,ny,ds]=EleBasedGL(CtrlVar,MUA,GF,DoPlots)
%
% Outputs:
%
%  xGLele, yGLele : ordered list of all corner points of GL edges, with NaN.
%                   between individual GL segments. Useful for plotting GL.
% xEdge, yEdge    : x y coordinates of centre points of each GL edge.
% nx, ny          : unit normal to the GL edges.
% ds              : length of each GL edge.
% FB              : end nodes of each edge.
% triGR           : triangulation of the grounded elements (matlab triangulation
%                   object)
%
% To plot FE mesh of grounded elements, the grounding line, and normals to the
% grounding line:
%
% Example:
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','CtrlVar','GF')
%   [xGLele,yGLele,triGR,FB,xEdge,yEdge,nx,ny,ds]=EleBasedGL(CtrlVar,MUA,GF) ; 
%
%   figure ;
%   triplot(triGR)                            % Plot the triangulation of grounded elements
%   hold on
%   quiver(xEdge,yEdge,nx,ny) ; axis equal ;  % Plot the normals to each grounding-line edges
%   plot(xGLele,yGLele,'r')                   % Plot the grounding line edges
%   plot(xEdge,yEdge,'.g')                    % Plot centre points of each grounding-line edge
%
%
%  
% See also PlotGroundingLines,   GLgeometry
%
%%


if isempty(CtrlVar)
    CtrlVar.GLthreshold=0.5;
end

if nargin<4 || isempty(DoPlots)
    DoPlots=0;
end

GF.ele=Nodes2EleMean(MUA.connectivity,GF.node);

if isfield(CtrlVar,'GLthreshold')
    I=GF.ele>CtrlVar.GLthreshold;
else
    I=GF.ele>0.5 ;
end

triGR=CreateFEmeshTriRep(MUA.connectivity(I,:),MUA.coordinates);
FB=freeBoundary(triGR);

xa=MUA.coordinates(FB(:,1),1);  ya=MUA.coordinates(FB(:,1),2);
xb=MUA.coordinates(FB(:,2),1);  yb=MUA.coordinates(FB(:,2),2);

ds=sqrt((xa-xb).^2+(ya-yb).^2);

nx=ya-yb ; ny=xb-xa ; temp=sqrt(nx.*nx+ny.*ny); nx=nx./temp ; ny=ny./temp;
xEdge=(xa+xb)/2; yEdge=(ya+yb)/2;

if ~(isfield(GF,'xGLele') && isfield(GF,'yGLele') && ~isempty(GF.xGLele) && ~isempty(GF.yGLele))
    [xGLele,yGLele]=LineUpEdges2(CtrlVar,xa,xb,ya,yb);
else
    xGLele=GF.xGLele ;
    yGLele=GF.yGLele;
end



if DoPlots
    figure ;
    triplot(triGR)
    hold on
    quiver(xEdge,yEdge,nx,ny) ; axis equal
    plot(xGLele,yGLele,'r')
    plot(xEdge,yEdge,'.g')
end

end
