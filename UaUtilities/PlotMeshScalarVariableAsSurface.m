function PatchObject=PlotMeshScalarVariableAsSurface(CtrlVar,MUA,Variable,AspectRatio,Col)

%%
%
% Plot a scalar nodal variable as a 3D surface.
%
%   PatchObject=PlotMeshScalarVariableAsSurface(CtrlVar,MUA,Variable,AspectRatio,Col)
%
%
% Note: The triangulation used when plotting will only be identical to the FE triangulation for
% 3-node elements.
%
% Example:
%
%    load('PIG-TWG-RestartFile.mat','CtrlVarInRestartFile','MUA','F')
%    Fig=FindOrCreateFigure("PIG-TWG Surface") ;
%    CtrlVar=CtrlVarInRestartFile ;
%    AspectRatio=50 ; 
%    PlotMeshScalarVariableAsSurface(CtrlVarInRestartFile,MUA,F.s,AspectRatio) ;
%    colormap(othercolor('Blues4',1024));
%
%
%
%  To change the colormap use for example:   
%
%   colormap(parula)  
%   colormap(othercolor('BuOr_12',1024));
%
%
% An example of how to specify  the colormap
%
%   Col=copper(numel(Variable));      
%   ColorIndex=Variable2ColorIndex(Variable);
%   Col(:,:)=Col(ColorIndex,:);
%
% And then give this as the input variable Col in the call.
%%



if nargin<5
    Col=[];
end

if nargin<4 || isempty(AspectRatio)
    AspectRatio=1;
end

TRI=TriFE(MUA.connectivity);
x=MUA.coordinates(:,1) ;
y=MUA.coordinates(:,2) ;

if isempty(Col)
    PatchObject=trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,Variable) ;
else
    PatchObject=trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,Variable,'FaceVertexCData',Col,'EdgeColor','none') ;
end


lighting phong ;

xlabel(CtrlVar.PlotsXaxisLabel) ;
ylabel(CtrlVar.PlotsYaxisLabel) ;
%zlabel('z (m a.s.l.)')
%colorbar ; title(colorbar,'(m)')
hold on



axis equal ; tt=daspect ;
daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)*CtrlVar.PlotXYscale/AspectRatio]);
axis tight



end