function [TRI,DT,LightHandle]=Plot_sbB(CtrlVar,MUA,s,b,B,TRI,DT,AspectRatio,ViewAndLight,LightHandle,sCol,bCol,BCol)

%%  Creates a perspective plot of s,b and B
%
% [TRI,DT,LightHandle]=Plot_sbB(CtrlVar,MUA,s,b,B,TRI,DT,AspectRatio,ViewAndLight,LightHandle)
%
%  TRI and DT are optional, can be empty.
%  If TRI and DT are available as outputs from a previous call, then give those
%  as input to following calls to speed things up.
%
%  Note: AspectRatio is not the actual aspect ratio between xy and z,
%        just a number that affects the aspect ratio.
%        To see what the aspect ratio is use the matlab commant `daspect'
%
% ViewAndLight(1)=AZ
% ViewAndLight(2)=EL
%
%
% CtrlVar.ThicknessCutOffForPlotting  :  Ice only plotted as ice if thickness greater than this.
%
% Examples:
%  
%  Plot_sbB(CtrlVar,MUA,s,b,B);
%
%  Plot_sbB(CtrlVar,MUA,Meas.s,Priors.Trueb,Priors.B,[],[],[],[],[],[1 0 0],[0 1 0],[0 0 1]) ;
%
%
% Note: TRI and DT are now calculated more efficiently and there is no longer
% any noticable gain in speed by giving those as inputs.
%%

x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ;

if nargin<10
    LightHandle=[];
end

if nargin<6 || isempty(TRI)
    %[TRI,DT]=CreateTRI(MUA);
    TRI=TriFE(MUA.connectivity);
end

if nargin<11
    sCol=[]; bCol=[] ; BCol=[];
end

if nargin<12
bCol=[];
end

if nargin<13
    BCol=[];
end


if nargin<8  || isempty(AspectRatio)
    AspectRatio=1;
end

if nargin<9 || isempty(ViewAndLight)
    ViewAndLight(1)=-70 ;  ViewAndLight(2)=20 ;
    ViewAndLight(3)=-45 ;  ViewAndLight(4)=50;
end

hold off


if ~isempty(s) && ~isempty(b)
    h=s-b;
else
    h=[] ;
end


if isfield(CtrlVar,'ThicknessCutOffForPlotting')
    I=h>CtrlVar.ThicknessCutOffForPlotting;
else
    I=h>2*CtrlVar.ThickMin;
end

if ~isempty(s)
    if isempty(sCol)
        sCol=copper(numel(s));
        ColorIndex=Variable2ColorIndex(s);
        sCol(:,:)=sCol(ColorIndex,:);
        sCol(I,:)=zeros(numel(find(I)),3)+1;
    end
end

if ~isempty(b)
    if isempty(bCol)
        bCol=copper(numel(s));
        ColorIndex=Variable2ColorIndex(b); bCol(:,:)=bCol(ColorIndex,:);
        bCol(I,:)=zeros(numel(find(I)),3)+1;
    end
end

if ~isempty(B)
    if isempty(BCol)
        BCol=copper(numel(s));
        ColorIndex=Variable2ColorIndex(B); BCol(:,:)=BCol(ColorIndex,:);
        BCol(I,:)=zeros(numel(find(I)),3)+1;
    end
end




if ~isempty(s)
    trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,s,'FaceVertexCData',sCol,'EdgeColor','none') ;

end

hold on
if ~isempty(b)
    trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,b,'FaceVertexCData',bCol,'EdgeColor','none') ;
end

if ~isempty(B)
    trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,B,'FaceVertexCData',BCol,'EdgeColor','none') ;
end



if ishandle(LightHandle)
    LightHandle=lightangle(LightHandle,ViewAndLight(3),ViewAndLight(4)) ;
else
    LightHandle=lightangle(ViewAndLight(3),ViewAndLight(4)) ;
end

lighting phong ;

xlabel(CtrlVar.PlotsXaxisLabel) ;
ylabel(CtrlVar.PlotsYaxisLabel) ;
zlabel('z (m a.s.l.)')
%colorbar ; title(colorbar,'(m)')
hold on

title(sprintf('t=%f (yr)',CtrlVar.time))
axis equal ; tt=daspect ;
daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)*CtrlVar.PlotXYscale/10/AspectRatio]);
axis tight
hold off

view(ViewAndLight(1),ViewAndLight(2));

end