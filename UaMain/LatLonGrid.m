



function LatLonGrid(X,Y,lat,lon,options)


%%
%
% plots lat lon grid on top of an existing figure
%
% X, Y          : x,y locations provided as meshgrid or ndgrid fields
%
% lat , lon     : grid of lat and lon values for the X,Y locations
%
% Note:
%
% To calculate lat lon from X,Y use something like:
%
%   [lat,lon]=psn2ll(X,Y);
%
%
% Example:
%
%   x=linspace(min(F.x),max(F.x),100) ; y=linspace(min(F.y),max(F.y),100)  ; 
%   [X,Y]=ndgrid(x,y) ; 
%   [lat,lon]=psn2ll(X,Y); 
%   hold on ; LatLonGrid(X/1000,Y/1000,lat,lon,LineColor=[0.5 0.5 0.5],LabelSpacing=200,LevelStepLat=5,LevelStepLon=10);
%
%%


arguments

    X    (:,:) double
    Y    (:,:) double
    lat  (:,:) double
    lon  (:,:) double

    options.LineColor double = [0.5 , 0.5 , 0.5]
    options.LabelSpacing double = 200
    options.FontSize double = 9
    options.LineStyle  string = "--"
    options.LevelStepLat = [];
    options.LevelStepLon = [];

end

% Guessing a sensible step between contour lines, if the user does not specify this in the call
if  isempty(options.LevelStepLat )
    latRange=max(lat(:))-min(lat(:));
    if latRange <0.5
        options.LevelStepLat = 0.1;
    elseif latRange<1.0
        options.LevelStepLat = 0.2;
    elseif latRange<10.0
        options.LevelStepLat = 1;
    else
        options.LevelStepLat = 5;
    end
end

if  isempty(options.LevelStepLon )
    lonRange=max(lon(:))-min(lon(:));
    if lonRange<0.5
        options.LevelStepLon = 0.1;
    elseif lonRange<1.0
        options.LevelStepLon = 0.2;
    elseif lonRange<10.0
        options.LevelStepLon = 1;
    else
        options.LevelStepLon = 5.0;
    end
end

hold on
[clat,hlat]=contour(X,Y,lat,LabelFormat=@mylabelfunLat,LineColor=options.LineColor,linestyle=options.LineStyle,LevelStep=options.LevelStepLat);
set(hlat,'ShowText','on','TextStep',get(hlat,'LevelStep')*2,'LabelSpacing',options.LabelSpacing)

clabel(clat,hlat,fontsize=options.FontSize) ;

[clon,hlon]=contour(X,Y,lon,LabelFormat=@mylabelfunLon,LineColor=options.LineColor,linestyle=options.LineStyle,LevelStep=options.LevelStepLon);
set(hlon,'ShowText','on','TextStep',get(hlon,'LevelStep')*2,'LabelSpacing',options.LabelSpacing)

clabel(clon,hlon,fontsize=options.FontSize)


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