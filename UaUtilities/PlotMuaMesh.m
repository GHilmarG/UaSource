function PlotMuaMesh(CtrlVar,MUA,ElementList,varargin)

%%
%
%   PlotMuaMesh(CtrlVar,MUA,ElementList,varargin)
%
% The only essential input is MUA, the others are optional.
%
% varargin is passed onto PlotFEmsh and then onto triplot.
%
% *Examples:*
%
% Plot Mesh:
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','CtrlVar')
%   figure ; PlotMuaMesh([],MUA)
%
% Plot the first 10000 elements in black:
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','CtrlVar')
%   figure ; PlotMuaMesh(CtrlVar,MUA,1:10000)
%
% Plot every 10th element in black and the nodes of those elements in black and do not plot the MeshBoundaryCoordinates
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','CtrlVar')
%   CtrlVar.PlotNodes=1; CtrlVar.NodeColor='r';
%   CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
%   figure ; PlotMuaMesh(CtrlVar,MUA,1:10:size(MUA.connectivity,1))
%
%
% Plot elements 3 to 10 in green and their nodes in red.
% Label both nodes and elements. Do not plot the MeshBoundaryCoordinates
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','CtrlVar')
%   CtrlVar.PlotNodes=1; CtrlVar.NodeColor='r';
%   CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
%   CtrlVar.PlotEleLabels=1;
%   CtrlVar.PlotNodalLabels=1;
%   CtrlVar.PlotNodesSymbol='*';
%   CtrlVar.PlotNodesSymbolSize=10;
%   figure ; PlotMuaMesh(CtrlVar,MUA,3:10,'g')
%
%%


if isempty(CtrlVar)
    CtrlVar.PlotLabels=0;
    CtrlVar.MeshColor='k';
    CtrlVar.NodeColor='k';
    CtrlVar.PlotXYscale=1;
    CtrlVar.PlotNodesSymbolSize=3;
    CtrlVar.PlotNodesSymbol='o';
    CtrlVar.PlotNodes=1;
    CtrlVar.time=NaN;
    CtrlVar.FEmeshPlotTitle=[];
    CtrlVar.PlotFEmeshAndSaveMesh=0;
    CtrlVar.PlotsXaxisLabel='x';
    CtrlVar.PlotsYaxisLabel='y';
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
end




if ~isfield(CtrlVar,'WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo')
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
end

if nargin<3  || isempty(ElementList)
    
    ElementList=1:MUA.Nele;
end

PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar,ElementList,varargin{:})



if CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo
    hold on
    PlotGmshGeometryDefinition(CtrlVar);
end

end