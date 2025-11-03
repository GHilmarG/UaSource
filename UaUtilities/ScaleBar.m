



function [ScaleLineHandle,ScaleTextHandle]=ScaleBar(Scale,ScaleLabel,ScalePosition,ScaleWidth,ScaleColor)

%%
%
% Creates a (very simple) scale bar to put on maps.
%
%
% Example:
%
%   ScaleBar();     % uses defaults
%
% 
% Create a 100 km scale-bar at relative x,y location of 0.8,0.8, i.e. towards the upper right corner of the figure
%
%   ScaleBar(100,"100 km",[0.8 0.8);
%
% See also:
%
%   PlotLatLonGrid
%%

xlim = get(gca,'xlim');
ylim = get(gca,'ylim');

if nargin==0

    if diff(xlim)<200

        Scale=10;
        ScaleLabel="10 km";

    elseif diff(xlim)<350

        Scale=50;
        ScaleLabel="50 km";

    elseif diff(xlim)< 1000

        Scale=100;
        ScaleLabel="100 km";

    elseif diff(xlim)< 2000

        Scale=500;
        ScaleLabel="500 km";


    else
        Scale=1000;
        ScaleLabel="1000 km";


    end


    ScalePosition(1)=0.025;
    ScalePosition(2)=0.05;
    ScaleWidth=3;
    ScaleColor="k";

end

x1=xlim(1)+(xlim(2)-xlim(1))*ScalePosition(1);
x2=x1+Scale;

y1=ylim(1)+(ylim(2)-ylim(1))*ScalePosition(2);


hold on
ScaleLineHandle=line([x1 x2],[y1 y1],9999*[1 1],color=ScaleColor,linewidth=ScaleWidth);

ScaleTextHandle=text(mean([x1 x2]),y1,ScaleLabel,...
    horizontalalignment="center",...
    verticalalignment="bottom",FontSize=12);

if nargout==0
    clear ScaleLineHandle
end


end