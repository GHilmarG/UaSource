function [Lat,Lon,X0,Y0,Clat,hlat,Clon,hlon]=PlotLatLonGrid(scale,dlat,dlon,LabelSpacing,Colour,isCircumpolar)

% Plots a lat lon grid
%
% [Lat,Lon,X0,Y0,Clat,hlat,Clon,hlon]=PlotLatLonGrid(scale,dlat,dlon,LabelSpacing)
%
%
%  hold on ; PlotLatLonGrid(1000,5,30,2000,1);
%


tt=axis;
xmin=tt(1) ; xmax=tt(2) ; ymin=tt(3) ; ymax=tt(4) ;

if nargin<2 || isempty(dlat)
    dlat=5;
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


if nargin<6
    isCircumpolar=0;
end

lcol='k';



[X0,Y0]=meshgrid(linspace(xmin,xmax,100),linspace(ymin,ymax,100));

[Lat,Lon]=pol_to_geog_wgs84_71S(X0*scale,Y0*scale);

if isCircumpolar
    I=Lat>-62; Lon(I)=NaN ; Lat(I)=NaN;
    I=Lat>-64.9;  Lon(I)=NaN;
    I=Lat<-85.1 ; Lon(I)=NaN;
    I=Lat<-86 ; Lat(I)=NaN ; 
    I=Lon<-175 ; Lon(I)=Lon(I)+360;
    I=Lon<-170 ; Lon(I)=NaN;
end


hold on
[Clat,hlat]=contour(X0,Y0,Lat,[-90:dlat:0],'LineColor',lcol);
set(hlat,'ShowText','on','TextStep',get(hlat,'LevelStep')*2,'LabelSpacing',LabelSpacing)

[Clon,hlon]=contour(X0,Y0,Lon,[-180+dlon:dlon:180],'LineColor',lcol);
set(hlon,'ShowText','on','TextStep',get(hlon,'LevelStep')*2,'LabelSpacing',LabelSpacing)

hlon.LineColor=Colour ; 
hlat.LineColor=Colour ; 
clabel(Clat,hlat,'Color',Colour)
clabel(Clon,hlon,'Color',Colour)




end
