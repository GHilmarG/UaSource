function [Lat,Lon,X0,Y0,Clat,hlat,Clon,hlon]=PlotLatLonGrid(scale,dlat,dlon,LabelSpacing,Colour)

%% 
% Plots a lat lon grid
% [Lat,Lon,X0,Y0,Clat,hlat,Clon,hlon]=PlotLatLonGrid(scale,dlat,dlon,LabelSpacing,Colour)
%
%      scale :  distance units, if using m set to 1 (default), if using km set to 1000
% dlat, dlon :  spacing between lat and lon lines (default is dlat=2, dlon=5)
%
% on return
%
% Lat, Lon the meshgrid with lat lon values
% Clat and Clon are the contour matrixes
% hlat and hlon the contour objects
%
% Examples:
% 
% PlotLatLonGrid(1000)   ; 
% xy axis in the units of km
%
% PlotLatLonGrid(1000,0.5,1,[],'r')   ; 
% xy axis in the units of km and lat lon, lines in red
%

tt=axis;
xmin=tt(1) ; xmax=tt(2) ; ymin=tt(3) ; ymax=tt(4) ;

if nargin<1 || isempty(scale)
    scale=1;
end

if nargin<2 || isempty(dlat)
    dlat=2;
end

if nargin<3 || isempty(dlon)
    dlon=5;
end

if nargin<4 || isempty(LabelSpacing)
    LabelSpacing=400;
end

if nargin<4 || isempty(Colour)
   Colour='black'; 
end

lcol='k';



[X0,Y0]=meshgrid(linspace(xmin,xmax,200),linspace(ymin,ymax,200));

[Lat,Lon]=pol_to_geog_wgs84_71S(X0*scale,Y0*scale);


hold on
[Clat,hlat]=contour(X0,Y0,Lat,[-90:dlat:0],'LineColor',lcol);
set(hlat,'ShowText','on','TextStep',get(hlat,'LevelStep')*2,'LabelSpacing',LabelSpacing)

%[dx,dy]=gradient(Lon);  Lon(abs(dx)>100)=180 ; 
%I=abs(X0)<100 & Y0<0;
%Lon(I)=abs(Lon(I));

%Lon=abs(Lon);

Lon(Lat<-85)=NaN;

I=X0<0 & Y0<0 & Lon<-175;

Lon(I)=NaN;

levels=-180+dlon:dlon:180;
[Clon,hlon]=contour(X0,Y0,Lon,levels,'LineColor',lcol);
set(hlon,'ShowText','on','TextStep',get(hlon,'LevelStep')*2,'LabelSpacing',LabelSpacing)

hlon.LineColor=Colour ; 
hlat.LineColor=Colour ; 
clabel(Clat,hlat,'Color',Colour)
clabel(Clon,hlon,'Color',Colour)


% [Clon,hlon]=contour(X0,Y0,abs(Lon),levels,'LineColor',lcol);
% set(hlon,'ShowText','on','TextStep',get(hlon,'LevelStep')*2,'LabelSpacing',LabelSpacing)
% 




end
