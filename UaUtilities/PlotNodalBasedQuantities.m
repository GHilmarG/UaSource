function [PatchObject,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,coordinates,Value,CtrlVar,varargin)

%% Plots nodal quantities as a patch
% [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,coordinates,Value,CtrlVar,varargin)
% plots nodal-based quantities in a map-plane view
%
% uses patch, varargin is passed on to patch
%
% tri is a 3-node triangulation. if tri is given as connectivity for 6 or 10 node elemets
% then a 3-node triangulation is created prior to plotting, and this triangulation is returned as tri
%
% Examples:
% Plot surface (s):
% PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,s,CtrlVar)
%
% Plot surface (s) and then bed (b). Use the output `tri' from the first call
% as an input to the second call.
%  [~,~,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,s,CtrlVar)
%  [~,~,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,b,CtrlVar)
%
%

%% Check inputs
if nargin<4 || isempty(CtrlVar)
    CtrlVar.PlotXYscale=1;
    CtrlVar.PlotsXaxisLabel=' ';
    CtrlVar.PlotsYaxisLabel=' ';
end

if ~isfield(CtrlVar,'PlotXYscale') 
    CtrlVar.PlotXYscale=1;
end

if ~isfield(CtrlVar,'PlotsXaxisLabel') 
    CtrlVar.PlotsXaxisLabel=' ';
end

if ~isfield(CtrlVar,'PlotsYaxisLabel') 
    CtrlVar.PlotsYaxisLabel=' ';
end

%%


[Nele,nod]=size(tri);


if nod~=3 
    tri=TriFE(tri);
end


Nvalues=length(Value);
Nnodes=size(coordinates,1);

if Nnodes~=Nvalues
    error('Number of nodal values to plot must be equal to number of nodes.')
end


%FigHandle=trisurf(tri,coordinates(:,1)/CtrlVar.PlotXYscale,coordinates(:,2)/CtrlVar.PlotXYscale,Value,'EdgeColor','none') ;
%view(0,90)
%tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);

PatchObject=patch('faces',tri,'vertices',coordinates/CtrlVar.PlotXYscale,...
    'FaceVertexCData',Value,'CDataMapping','scaled','EdgeColor','none','FaceColor','interp',varargin{:}) ;

axis equal
ColorbarHandel=colorbar;
xlabel(CtrlVar.PlotsXaxisLabel)  ; ylabel(CtrlVar.PlotsYaxisLabel) ;


return
