function [Lat,Lon,X0,Y0,Clat,hlat,Clon,hlon,ax1,ax2]=PlotLatLonGrid(scale,dlat,dlon,LabelSpacing,Colour,isCircumpolar)

%%
% Plots a lat lon grid 
%
% This is written for the Antarctica setting using polar stereographic coordinate system.
%
% 
%
% [Lat,Lon,X0,Y0,Clat,hlat,Clon,hlon]=PlotLatLonGrid(scale,dlat,dlon,LabelSpacing,Colour,isCircumpolar)
%
%
% NOTE: Do not use figure zoom after this or the lat/lon lin will get misaligned!
%       Despite best atempts I have not been able to link the axis and get the right behaviour.
%
% Example:
%
%   load('PIG-TWG-RestartFile.mat','CtrlVarInRestartFile','MUA','F')
%   Fig=FindOrCreateFigure("PIG-TWG lat/lon") ;
%   CtrlVar=CtrlVarInRestartFile ;
%   cbar=UaPlots(CtrlVar,MUA,F,"-speed-") ; 
%   hold on ; 
%   [~,~,~,~,~,hlat,~,hlon]=PlotLatLonGrid(1000)   ; % often the colormap will have to be redefined after this call
%   axis([-2000 -1000 -900 100])
%   hlat.LineStyle="--"; hlon.LineStyle="--";
%   clim([0 4000])
%   ModifyColormap;
%
%
%

%%

fig = gcf;
ax1 = fig.CurrentAxes ; 
tt=axis;

xmin=tt(1) ; xmax=tt(2) ; ymin=tt(3) ; ymax=tt(4) ;


%% create new axes for the lat/lon lines  (never got this to work)
%     ax2=axes ;
%     
%     ax2.Visible = 'off';
%     ax2.XTick = [];
%     ax2.YTick = [];
%     hold on
%     ax2.Position=ax1.Position;
%     ax2.XLim=ax1.XLim;
%     ax2.YLim=ax1.YLim;
ax2=[] ; 
%%

if nargin<2 || isempty(dlat)
    dlat=5;
end

if nargin<3 || isempty(dlon)
    dlon=10;
end

if nargin<4 || isempty(LabelSpacing)
    LabelSpacing=400;
end

if nargin<4 || isempty(Colour)
   Colour='black'; 
end


if nargin<6
    isCircumpolar=0;
end

lcol='k';



[X0,Y0]=meshgrid(linspace(xmin,xmax,400),linspace(ymin,ymax,400));

[Lat,Lon]=pol_to_geog_wgs84_71S(X0*scale,Y0*scale);

if isCircumpolar
    I=Lat>-62; Lon(I)=NaN ; Lat(I)=NaN;
    I=Lat>-64.9;  Lon(I)=NaN;
    I=Lat<-85.1 ; Lon(I)=NaN;
    I=Lat<-86 ; Lat(I)=NaN ; 
    I=Lon<-171 ; Lon(I)=Lon(I)+360;
    I=Lon<-170 ; Lon(I)=NaN;
end


hold on
[Clat,hlat]=contour(ax1,X0,Y0,-Lat,[0:dlat:90],'LineColor',lcol,"LabelFormat","%2.0fS");
set(hlat,'ShowText','on','TextStep',get(hlat,'LevelStep')*2,'LabelSpacing',LabelSpacing)


[Clon,hlon]=contour(ax1,X0,Y0,Lon,[-180+dlon:dlon:185],LineColor=lcol,LabelFormat=@mylabelfun);
set(hlon,'ShowText','on','TextStep',get(hlon,'LevelStep')*2,'LabelSpacing',LabelSpacing)


hlon.LineColor=Colour ; 
hlat.LineColor=Colour ; 
clabel(Clat,hlat,Color=Colour,fontsize=9);
clabel(Clon,hlon,Color=Colour,fontsize=9)


%linkaxes([ax1,ax2],"xy") ; %  For some reason this is not having the desired effect...?!
%fig.CurrentAxes = ax1;
%ax2.Position=ax1.Position;
% revert back to original axes




    function labels=mylabelfun(vals)


    labels= vals +"E" ; 
    I=vals<0  ;  labels(I) = -vals(I) +"W" ; 
    I=vals==0 ;  labels(I) = vals(I)  ; 
    

    end







end
