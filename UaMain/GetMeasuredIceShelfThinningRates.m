function [Fdhdt,dhdt]=GetMeasuredIceShelfThinningRates(xps,yps,dhdzGL)

%%

if nargin<3
    dhdzGL=NaN;
end

AntarcticGlobalDataSets=getenv('AntarcticGlobalDataSets');


locdir=pwd;


cd(AntarcticGlobalDataSets);
cd IceShelfThinningRates/dzdt_v3/dzdt_v3/
%%


data=load('ice_shelf_dzdt_v3.mat');
[x,y]=ll2xy(data.lat,data.lon);
x=x(:) ; y=y(:) ; dzdt=data.dzdt_raw(:);


% cd(AntarcticGlobalDataSets);
% cd IceShelfThinningRates/dzdt/dzdt/
% load('ice_shelf_dzdt_v01.mat','lat','lon','dzdt')
% 
% [x,y]=ll2xy(lat,lon);
% x=x(:) ; y=y(:) ; dzdt=dzdt(:);

cd(locdir)

% 
% Interpolation='ReplaceNaNWithInterpolatedValues';
% 
% switch Interpolation
% 
%         
%     case 'ReplaceNaNWithInterpolatedValues'
%         
%         
%         I=~isnan(dzdt); x=x(I) ; y=y(I); dzdt=dzdt(I); 
%         
%         % Excluding NaN. This has the effect that the interpolation/extrapolation gives values over the whole area of each and every ice shelf
%         % If I keep the NaNs, then I get NaN in interpolated values that I can then later replace with zeros
%         
% end

if ~isnan(dhdzGL)
    [xGLbob,yGLbob]=ReadBindschadlerGroundingLine;
    R=sqrt(x.*x+y.*y);
    Rgl=sqrt(xGLbob.*xGLbob+yGLbob.*yGLbob);
    I=Rgl<min(R);
    
    % add fictious data points along grounding line where no data is available
    % set these points to some values representing possible range
    
    
    xDataAdded=xGLbob(I) ;  yDataAdded=yGLbob(I) ;
    DataAdded=yDataAdded*0+dhdzGL;
    
    
    x=[x(:) ; xDataAdded];
    y=[y(:) ; yDataAdded];
    dzdt=[dzdt(:) ; DataAdded];
    I=~isnan(dzdt); x=x(I) ; y=y(I); dzdt=dzdt(I);
    
end

Fdhdt=scatteredInterpolant(x,y,dzdt,'natural','nearest');

if nargin>=2
    dhdt=Fdhdt(xps,yps);
else
    dhdt=[];
end


end