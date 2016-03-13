function GmshPointsInLoop=PlotGmshGeometryDefinition(CtrlVar)

% Plots the geometry (topology) of the mesh as given as input to gmsh
% when using input format 2 (CtrlVar.GmshInputFormat=2)
%
% GmshPointsInLoop=PlotGmshGeometryDefinition(CtrlVar)
%
%
%

if CtrlVar.GmshInputFormat==2
    
    GmshPointsInLoop=FindGmshPointsInLoop(CtrlVar);
    
    
    for I=1:numel(GmshPointsInLoop)
        plot(GmshPointsInLoop{I}(:,1)/CtrlVar.PlotXYscale,GmshPointsInLoop{I}(:,2)/CtrlVar.PlotXYscale,'LineWidth',2)
    end
    
else
    
    fprintf('PlotGmshGeometryDefinition: Not using input format 2 to gmsh')
    
end


end
