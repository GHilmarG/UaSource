function PatchObject=PlotMeshScalarVariableAsSurface(CtrlVar,MUA,Variable,AspectRatio,Col)


%
%
%  To change the colormap use for example:   
%
%   colormap(parula)  
%   colormap(othercolor('BuOr_12',1024));
%
%

if nargin<5
    Col=[];
end

if nargin<4 || isempty(AspectRatio)
    AspectRatio=1;
end

TRI=TriFE(MUA.connectivity);
x=MUA.coordinates(:,1) ;
y=MUA.coordinates(:,2) ;


% An example of how to change the colormap
% if ~isempty(Variable)
%     if isempty(Col)
%         Col=copper(numel(Variable));
%         ColorIndex=Variable2ColorIndex(Variable);
%         Col(:,:)=Col(ColorIndex,:);
%     end
% end


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