


function [lat,lon]=xy2LatLonIceland(x,y)

mstruct=defaultm('lambertstd');
mstruct.geoid = almanac('earth','geoid','m','grs80');
%mstruct.geoid = referenceEllipsoid('wgs84','meter');
mstruct.mapparallels=[64.25 65.75];
mstruct.origin=[65.0 -19.0];
mstruct.falsenorthing=500000;
mstruct.falseeasting=500000;
mstruct=defaultm(mstruct);

%[X,Y]=mfwdtran(mstruct,Lat,Lon);
%[X,Y]=mfwdtran(mstruct,Lat,Lon);


[lat,lon]=projinv(mstruct,x,y) ;



end
