

function [InfluxNodes,OutfluxNodes,InOutBoundary]=InfluxOutfluxNodes(CtrlVar,MUA,F)

%%
%
% Finds nodes along mesh boundary where ice flow is directed either into
% the boundary elements, or away.
%
% InfluxNodes       : list of nodes where the velocity vector points inwards
% OutfluxNodes      : list of nodes where the velocity vector points outwards
% 
% InOutBoundary     : a logical list of in/out nodes within MUA.Boundary.Nodes
%
%   InfluxNodes=MUA.Boundary.Nodes(InOutBoundary);
%   OutfluxNodes=MUA.Boundary.Nodes(~InOutBoundary);
%
%%



[nx,ny,xn,yn,Nx,Ny] = CalcEdgeAndNodalNormals(MUA.connectivity,MUA.coordinates,MUA.Boundary.Edges);
P=F.ub(MUA.Boundary.Nodes).*Nx(MUA.Boundary.Nodes)+F.vb(MUA.Boundary.Nodes).*Ny(MUA.Boundary.Nodes);
InOutBoundary=P<0 ;

InfluxNodes=MUA.Boundary.Nodes(InOutBoundary);
OutfluxNodes=MUA.Boundary.Nodes(~InOutBoundary);

return

%%  Example of how to plot some of the outputs

FindOrCreateFigure("Influx Nodes")
hold off
PlotMuaBoundary(CtrlVar,MUA);
hold on ; 
plot(F.x(MUA.Boundary.Nodes(:,1)) /CtrlVar.PlotXYscale,F.y(MUA.Boundary.Nodes(:,1))/CtrlVar.PlotXYscale,"*b")
hold on


% QuiverColorGHG(MUA.coordinates(MUA.Boundary.Nodes,1),MUA.coordinates(MUA.Boundary.Nodes,2),...
%     Nx(MUA.Boundary.Nodes),Ny(MUA.Boundary.Nodes),CtrlVar);

QuiverColorGHG(MUA.coordinates(MUA.Boundary.Nodes,1),MUA.coordinates(MUA.Boundary.Nodes,2),...
    F.ub(MUA.Boundary.Nodes),F.vb(MUA.Boundary.Nodes),CtrlVar);

hold on

plot(F.x(MUA.Boundary.Nodes(InOutBoundary))/CtrlVar.PlotXYscale,F.y(MUA.Boundary.Nodes(InOutBoundary))/CtrlVar.PlotXYscale,"or",MarkerFaceColor="r")

title("Influx nodes shown as red circles")

end


%%

