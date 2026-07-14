function [x,y]=LatLon2xyIceland(lat,lon)

%%
%
% Forward projection from (x,y) Icelandic geographic coordinate system to lat lon.
%
% Note: Many data sets seem to be provided in UTM zone ??N instead. In that case one can map to lat lon as follows:
% 
%  mstruct = projcrs(326??);  % e.g. EPSG 32627 is the standard EPSG code for WGS84 / UTM zone 27N
%  [Lat, Lon] = projinv(mstruct, X, Y);
%
% the Icelandic system is know as: EPSG 3057
%
% requires the MATLAB mapping toolbox
%
% see also: xy2LatLonIceland
%%

mstruct=defaultm('lambertstd');
mstruct.geoid = almanac('earth','geoid','m','grs80');
%mstruct.geoid = referenceEllipsoid('wgs84','meter');
mstruct.mapparallels=[64.25 65.75];
mstruct.origin=[65.0 -19.0];
mstruct.falsenorthing=500000;
mstruct.falseeasting=500000;
mstruct=defaultm(mstruct);




[x,y]=projfwd(mstruct,lat,lon) ;



end
