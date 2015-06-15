function [TRI,DT]=Plot_sbB(CtrlVar,MUA,s,b,B,TRI,DT,AspectRatio,ViewAndLight)

% Creates a perspective plot of s,b and B
%  TRI and DT are optional, can be empty.
%
%

persistent light_handle

x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ;


if nargin<6 || isempty(TRI)
    [TRI,DT]=CreateTRI(MUA);
end

if nargin<8  || isempty(AspectRatio)
    AspectRatio=1;
end

if nargin<9
    ViewAndLight(1)=-70 ;  ViewAndLight(2)=20 ;
    ViewAndLight(3)=-45 ;  ViewAndLight(4)=50;
end

hold off
if ~isempty(s)
    trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,s,'EdgeColor','none') ;
end

hold on
if ~isempty(b)
    trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,b,'EdgeColor','none') ;
end

if ~isempty(B)
    trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,B,'EdgeColor','none') ;
end


view(ViewAndLight(1),ViewAndLight(2))

if ishandle(light_handle)
    light_handle=lightangle(light_handle,ViewAndLight(3),ViewAndLight(4)) ;
else
    light_handle=lightangle(ViewAndLight(3),ViewAndLight(4)) ;
end

lighting phong ;

xlabel(CtrlVar.PlotsXaxisLabel) ;
ylabel(CtrlVar.PlotsYaxisLabel) ;

colorbar ; title(colorbar,'(m)')
hold on

title(sprintf('t=%#6.2f ',CtrlVar.time))
axis equal ; tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)*CtrlVar.PlotXYscale/10/AspectRatio]); 
axis tight
hold off

end