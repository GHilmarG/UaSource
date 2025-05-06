

function [InfluxNodes,OutfluxNodes,InOutBoundary]=InfluxOutfluxNodes(CtrlVar,MUA,F,options)

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

arguments
    CtrlVar struct
    MUA     struct
    F       {mustBeA(F,{'struct','UaFields','numeric'})}
    options.plot  logical = false
    options.afloat logical = false
    options.MinSpeed = 0 ;
end

[~,~,~,~,Nx,Ny] = CalcEdgeAndNodalNormals(MUA.connectivity,MUA.coordinates,MUA.Boundary.Edges);
P=F.ub(MUA.Boundary.Nodes).*Nx(MUA.Boundary.Nodes)+F.vb(MUA.Boundary.Nodes).*Ny(MUA.Boundary.Nodes);

InOutBoundary=P<0 ;

if options.afloat

   InOutBoundary = InOutBoundary &  F.GF.node(MUA.Boundary.Nodes)< 0.5 ; 

    
end


if options.MinSpeed > 0

  speed=sqrt(F.ub(MUA.Boundary.Nodes).^2+F.vb(MUA.Boundary.Nodes).^2) ; 

  InOutBoundary=InOutBoundary & speed > options.MinSpeed; 

end





InfluxNodes=MUA.Boundary.Nodes(InOutBoundary);
OutfluxNodes=MUA.Boundary.Nodes(~InOutBoundary);

if options.plot

    %%  Example of how to plot some of the outputs

    fIN=FindOrCreateFigure("Influx Nodes") ; clf(fIN);
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

    title("Selected influx nodes shown as red circles")

end

end


%%

