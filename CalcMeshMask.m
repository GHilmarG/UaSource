function Mask=CalcMeshMask(CtrlVar,MUA,NodalField,Threshold)

%%
%   Mask=CalcMeshMask(CtrlVar,MUA,NodalField,Threshold)
%
% 
% In general a (nodal) mask is an nodal array of values between 0 and 1 indicating
% if a the value of the NodalField at a given node is lower or higher than
% the specified Threshold.
% 
% Here the Mask is returned is defined slighly differently or as logical array where:
% 
% Mask.NodesOn  :  True for nodes belonging to elements where some of the nodal values are above and some below the Threshold
% Mask.NodesIn  :  True for nodes belonging to elements where all of the nodal values are above the Threshold
% Mask.NodesOut :  True for nodes belonging to elements where all of the nodal values are below the Threshold
% 
% A corresponding element mask is also returned. 
% 
% This definition of In/Out is sometimes also referred to as being 'strictly' above or
% 'strickly' below a given threshold value. 
%
%
% Example: Create a mask based on ice thickness above and below 100. 
%
%   Mask=CalcMeshMask(CtrlVar,MUA,F.h,100); 
% 
%   FindOrCreateFigure("Level Set Mask") ; 
%   CtrlVar.PlotNodes=0; CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
%   PlotMuaMesh(CtrlVar,MUA) ; 
%   hold on ; 
%   x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ; 
%   p1=plot(x(Mask.NodesIn)/CtrlVar.PlotXYscale,y(Mask.NodesIn)/CtrlVar.PlotXYscale,'or','DisplayName','In');
%   p2=plot(x(Mask.NodesOn)/CtrlVar.PlotXYscale,y(Mask.NodesOn)/CtrlVar.PlotXYscale,'og','DisplayName','On');
%   p3=plot(x(Mask.NodesOut)/CtrlVar.PlotXYscale,y(Mask.NodesOut)/CtrlVar.PlotXYscale,'ob','DisplayName','Out');
%   legend([p1 p2 p3])
% 
%%

narginchk(3,6)
nargoutchk(1,4)

if isempty(NodalField)
    return
end


EleValue=Nodes2EleMean(MUA.connectivity,NodalField);

CtrlVar.GLthreshold=Threshold ; % need this in GLgeometry
Mask.Level=CtrlVar.GLthreshold; 
[Mask.Geo,Mask.NodesOn,Mask.ElementsOn]=GLgeometry(MUA.connectivity,MUA.coordinates,NodalField,CtrlVar);

Mask.ElementsIn=   (EleValue>Mask.Level)  & ~Mask.ElementsOn;
Mask.ElementsOut= (EleValue<Mask.Level)  & ~Mask.ElementsOn;


Mask.NodesIn=false(MUA.Nnodes,1);
Mask.NodesOut=false(MUA.Nnodes,1);

Mask.NodesIn(MUA.connectivity(Mask.ElementsIn,:))=true;
Mask.NodesOut(MUA.connectivity(Mask.ElementsOut,:))=true;


Mask.NodesOut=Mask.NodesOut & ~Mask.NodesOn;
Mask.NodesIn=Mask.NodesIn & ~Mask.NodesOn;


return



end