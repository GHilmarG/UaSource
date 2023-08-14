function hTri=PlotMuaMesh(CtrlVar,MUA,ElementList,varargin)

%%
%
%   PlotMuaMesh(CtrlVar,MUA,ElementList,varargin)
%
% The only essential input is MUA, the other input variables are optional.
%
% varargin is passed onto PlotFEmsh and then onto triplot.
%
% hTri is a handle to the matlab triplot function
%
% *Examples:*
%
% Plot Mesh in red:
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','CtrlVar')
%   figure ; PlotMuaMesh([],MUA,[],'r')
%
% or
%
%   figure ; PlotMuaMesh(CtrlVar,MUA,[],color="r");
%
%
% Plot the first 10000 elements in black:
%
%   load('MUA-PIG-TWG-Example.mat','MUA','BCs','CtrlVar')
%   figure ; PlotMuaMesh(CtrlVar,MUA,1:10000)
%
% Plot every 10th element in black and the nodes of those elements in red and do not plot the MeshBoundaryCoordinates
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

if nargin<3  
    ElementList=1:MUA.Nele;
end



if ischar(ElementList) && nargin==3
    % silently ignore the fact that the user clearly did not read the comments and 
    % has given a character as third argument. Assume that he/she wanted
    % to specify some plotting options passed on to PlotFEmesh.
    varargin{1}=ElementList;
    ElementList=1:MUA.Nele;
end



hTri=PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar,ElementList,varargin{:}) ; 



if CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo
    hold on
    PlotGmshGeometryDefinition(CtrlVar);
end


[Emin,Emax,Emean,Emedian]=PrintInfoAboutElementsSizes(CtrlVar,MUA,print=false) ;


title(sprintf("#Ele=%i #Nodes=%i #nod=%i \n (max,mean,median,min)=(%g,%g,%g,%g) ",MUA.Nele,MUA.Nnodes,MUA.nod,Emax,Emean,Emedian,Emin))

if ~nargout   % A trick to suppress any function output if no output requested. No need to suppress output using ;
    clearvars hTri
end



end

