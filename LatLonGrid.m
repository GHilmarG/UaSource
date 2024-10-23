



function LatLonGrid(X,Y,lat,lon,options)


%%
%
% plots lat lon grid on top of an existing figure
%
% X, Y          : x,y locations provided as meshgrid or ndgrid fields
%
% lat , lon     : grid of lat and lon values for the X,Y locations
%
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

end



fig = gcf;


hold on 
[clat,hlat]=contour(X,Y,lat,LabelFormat=@mylabelfunLat,LineColor=options.LineColor,linestyle=options.LineStyle);
set(hlat,'ShowText','on','TextStep',get(hlat,'LevelStep')*2,'LabelSpacing',options.LabelSpacing)

clabel(clat,hlat,fontsize=options.FontSize) ;

[clon,hlon]=contour(X,Y,lon,LabelFormat=@mylabelfunLon,LineColor=options.LineColor,linestyle=options.LineStyle);
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