




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
% Inputes:
%
%  scale   ;    scales the x and the y axis. For example if the x,y units are meters, but you want to plot using km as a
%               distance units, set scale=1000
%
%  dlat, dlon : distance between lat and lon lines, in degrees
%
%  LableSpacing :   affects the spacing between labels, default is LableSpacing=400. Increase to increase spacing between lat
%                   lon labels.
%
% Colour:   color of the lat, lon lines
% 
% isCircumpolar:  set to true if the plot area is circumpolar, ie includes the pole itself.  
%
%
% Outputs:
%
%     Clat,hlat,Clon,hlon  : These are the contour matrices and the contour objects. The contour objects allow the properties
%     of the countour plots to be easily edited after the call.
%
%
%
% NOTE #1: Do not use figure zoom after this or the lat/lon lin will get misaligned!
%       Despite best atempts I have not been able to link the axis and get the right behaviour.
%
% NOTE #2: As of Matlab2023b, Note#1 is no longer of relevance, which is good news! 
%
%
% Example:
%
%   load('PIG-TWG-RestartFile.mat','CtrlVarInRestartFile','MUA','F')
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


%% Setting values in the case of no or limited input

% First try to figure out if this is circumpolar
if nargin< 6 || isempty(isCircumpolar)
    AxLimits=axis;
    isCircumpolar=AxLimits(1) < 0 && AxLimits(2) > 0  && AxLimits(3) <0 && AxLimits(4) > 0 ;
end

if isCircumpolar

    if nargin < 5
        Colour='black';
        if nargin< 4
            LabelSpacing=1000;
            if nargin< 3
                dlon=45;
                if nargin< 2
                    dlat=10;
                    if nargin==0
                        scale=1000 ;
                    end
                end
            end
        end
    end

else

    if nargin < 5
        Colour='black';
        if nargin< 4
            LabelSpacing=400;
            if nargin< 3
                dlon=10;
                if nargin< 2
                    dlat=2.5 ;
                    if nargin==0
                        scale=1000 ;
                    end
                end
            end
        end
    end

end

%%

% set some plausible values if user has not defined those already
if isCircumpolar && isempty(dlat) &&  isempty(dlon)  &&  isempty(LabelSpacing)

    dlat=10;
    dlon=45;
    LabelSpacing=200;

else


    if isempty(dlat)
        dlat=5;
    end

    if isempty(dlon)
        dlon=10;
    end

    if isempty(LabelSpacing)
        LabelSpacing=400;
    end

    if isempty(Colour)
        Colour='black';
    end


end

lcol='k';

climCopy=clim;


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


[Clat,hlat]=contour(ax1,X0,Y0,Lat,-90:dlat:90,LineColor=lcol,LabelFormat=@mylabelfunLat);


set(hlat,'ShowText','on','TextStep',get(hlat,'LevelStep')*2,'LabelSpacing',LabelSpacing)


[Clon,hlon]=contour(ax1,X0,Y0,Lon,-180+dlon:dlon:185,LineColor=lcol,LabelFormat=@mylabelfunLon);
set(hlon,'ShowText','on','TextStep',get(hlon,'LevelStep')*2,'LabelSpacing',LabelSpacing)


hlon.LineColor=Colour ;
hlat.LineColor=Colour ;
clabel(Clat,hlat,Color=Colour,fontsize=9);
clabel(Clon,hlon,Color=Colour,fontsize=9)


%linkaxes([ax1,ax2],"xy") ; %  For some reason this is not having the desired effect...?!
%fig.CurrentAxes = ax1;
%ax2.Position=ax1.Position;
% revert back to original axes


clim(climCopy) % set color axis limit to the value at the beginning of the call
               % this is done here because the contour functions above might change the existing limites

    function labels=mylabelfunLon(vals)

        % Degree=string(char(176))

        labels= vals +"°E" ;
        I=vals<0  ;  labels(I) = -vals(I) + "°W" ;
        I=vals==0 ;  labels(I) = vals(I)  ;


    end



    function labels=mylabelfunLat(vals)


        labels= vals +"°N" ;
        I=vals<0  ;  labels(I) = -vals(I) + "°S" ;
        I=vals==0 ;  labels(I) = vals(I)  ;


    end




end
