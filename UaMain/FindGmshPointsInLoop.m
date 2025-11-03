function GmshPointsInLoop=FindGmshPointsInLoop(CtrlVar)

% find the geometry (topology) of the mesh as given as input to gmsh
% when using input format 2 (CtrlVar.GmshInputFormat=2)
% 
% GmshPointsInLoop=PlotGmshGeometryDefinition(CtrlVar)
%
%
%

GmshPointsInLoop=[];

% plotting loops
if CtrlVar.GmshInputFormat==2
    
    for J=1:numel(CtrlVar.Gmsh.Loops)
        nod=[];
        temp={CtrlVar.Gmsh.Lines{abs(CtrlVar.Gmsh.Loops{J})}};
        for I=1:numel(temp)
            nod=[nod;temp{I}(:)];
        end
        %plot(CtrlVar.Gmsh.Points(nod,1)/CtrlVar.PlotXYscale,CtrlVar.Gmsh.Points(nod,2)/CtrlVar.PlotXYscale,'LineWidth',2)
        %hold on
        GmshPointsInLoop{J}=CtrlVar.Gmsh.Points(nod,:);
    end
    
else
    
    fprintf('PlotGmshGeometryDefinition: Not using input format 2 to gmsh')
    
end


end
