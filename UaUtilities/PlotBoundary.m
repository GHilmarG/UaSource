function PlotBoundary(Boundary,connectivity,coordinates,CtrlVar,varargin)

%% Plots/labels boundary nodes, edges and elements
%
% PlotBoundary(Boundary,connectivity,coordinates,CtrlVar,varargin)
%
% Examples:
% PlotBoundary(Boundary,connectivity,coordinates,CtrlVar)
% 
% PlotBoundary([],connectivity,coordinates,CtrlVar)
%
% if Boundary=[], it is created by a call to FindBoundary(connectivity,coordinates)
%
% For just plotting the boundary edges as, for example, a black line, do:
% figure ; PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,'k')
%
% Notes: 
% 
% Another option of just plotting the mesh boundaries is simply to do: 
% figure ; plot(MUA.coordinates(MUA.Boundary.Edges,1)/CtrlVar.PlotXYscale, MUA.coordinates(MUA.Boundary.Edges,2)/CtrlVar.PlotXYscale, 'k', 'LineWidth',2) ;
% in which case the `PlotBoundary' m-file is not used.
%
%
% See also: 
% PlotMuaBoundary
%%

if nargin < 4 || isempty(CtrlVar)
    CtrlVar.PlotLabels=0;
    CtrlVar.MeshColor='k';
    CtrlVar.NodeColor='k';
    CtrlVar.PlotXYscale=1;
    CtrlVar.PlotNodesSymbolSize=3;
    CtrlVar.PlotNodesSymbol='o';
    CtrlVar.PlotNodes=0;
    CtrlVar.PlotMesh=1;
    CtrlVar.time=NaN;
    CtrlVar.PlotBoundaryLabels=0;
    CtrlVar.PlotBoundaryElements=0;
    CtrlVar.PlotBoundaryNodes=0;
else
    if ~isfield(CtrlVar,'PlotLabels') ; CtrlVar.PlotLabels=0; end
    if ~isfield(CtrlVar,'MeshColor') ; CtrlVar.MeshColor='k'; end
    if ~isfield(CtrlVar,'NodeColor') ; CtrlVar.NodeColor='k'; end
    if ~isfield(CtrlVar,'PlotXYscale') ; CtrlVar.PlotXYscale=1; end
    if ~isfield(CtrlVar,'PlotNodesSymbol') ; CtrlVar.PlotNodesSymbol='o'; end
    if ~isfield(CtrlVar,'PlotNodesSymbolSize') ; CtrlVar.PlotNodesSymbolSize=3; end
    if ~isfield(CtrlVar,'PlotNodes') ; CtrlVar.PlotNodes=0; end
    if ~isfield(CtrlVar,'PlotMesh') ; CtrlVar.PlotMesh=1; end
    if ~isfield(CtrlVar,'time') ; CtrlVar.time=NaN; end
    if ~isfield(CtrlVar,'PlotBoundaryLabels'); CtrlVar.PlotBoundaryLabels=0;  end
    if ~isfield(CtrlVar,'PlotBoundaryElements'); CtrlVar.PlotBoundaryElements=0; end
    if ~isfield(CtrlVar,'PlotBoundaryNodes'); CtrlVar.PlotBoundaryNodes=0; end
end


if isempty(Boundary)
    [Boundary,~]=FindBoundary(connectivity,coordinates);
end


if ~isfield(Boundary,'x') || ~isfield(Boundary,'y')
    xa=coordinates(Boundary.Edges(:,1),1); xb=coordinates(Boundary.Edges(:,end),1);
    ya=coordinates(Boundary.Edges(:,1),2); yb=coordinates(Boundary.Edges(:,end),2);
    [Boundary.x,Boundary.y]=LineUpEdges2([],xa,xb,ya,yb);
    
end

plot(Boundary.x/CtrlVar.PlotXYscale,Boundary.y/CtrlVar.PlotXYscale,varargin{:}) ;

% plot boundary elements,  do not label nodes or elements here because
% the labelling will be incorrect
CtrlVar.PlotLabels=0;
if CtrlVar.PlotBoundaryElements
    PlotFEmesh(coordinates,connectivity(Boundary.FreeElements,:),CtrlVar)
end

hold on

if CtrlVar.PlotBoundaryNodes
    plot(coordinates(Boundary.EdgeCornerNodes,1)/CtrlVar.PlotXYscale,coordinates(Boundary.EdgeCornerNodes,2)/CtrlVar.PlotXYscale,'ok')
    plot(coordinates(Boundary.Nodes,1)/CtrlVar.PlotXYscale,coordinates(Boundary.Nodes,2)/CtrlVar.PlotXYscale,'xr')
end



% this gives correct labeling of nodes and elemetns
if CtrlVar.PlotBoundaryLabels
    labels = arrayfun(@(n) {sprintf(' N%d', n)}, Boundary.Nodes(:));
    text(coordinates(Boundary.Nodes,1)/CtrlVar.PlotXYscale,coordinates(Boundary.Nodes,2)/CtrlVar.PlotXYscale,labels)
    
    xEle=Nodes2EleMean(connectivity,coordinates(:,1));
    yEle=Nodes2EleMean(connectivity,coordinates(:,2));
    labels = arrayfun(@(n) {sprintf(' T%d', n)}, Boundary.FreeElements(:));
    text(xEle(Boundary.FreeElements)/CtrlVar.PlotXYscale,yEle(Boundary.FreeElements)/CtrlVar.PlotXYscale,labels,...
        'HorizontalAlignment', 'center', 'Color', 'blue');
    
end

%title(sprintf('boundary t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')
%axis equal tight

ax=gca; ax.DataAspectRatio=[1 1 1];


end