


function [lat,lon]=xy2LatLonIceland(x,y)

%%
% projects lat lon to a common xy system used in Iceland
%
%
% requires the MATLAB mapping toolbox
%
% 
% the Icelandic system is know as: EPSG 3057
% 
% 
% see also: LatLon2xyIceland
%
%
%
%%

mstruct=defaultm('lambertstd');
mstruct.geoid = almanac('earth','geoid','m','grs80');
%mstruct.geoid = referenceEllipsoid('wgs84','meter');
mstruct.mapparallels=[64.25 65.75];
mstruct.origin=[65.0 -19.0];
mstruct.falsenorthing=500000;
mstruct.falseeasting=500000;
mstruct=defaultm(mstruct);


[lat,lon]=projinv(mstruct,x,y) ;



end
