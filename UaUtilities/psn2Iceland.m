
function [x,y,u,v]=psn2Iceland(xpsn,ypsn,upsn,vpsn)

%% mapping from polar stereo-graphic north to the Icelandic coordinate system  EPSG 3057
%
%  maps the velocity field (upsn,vpsn) at locations (xpsn,ypsn) to the Icelandic system.
%
% The velocity rotation is a slight approximation and is calculated only at the centre of the Icelandic system, but this
% should be fine. 
%
% xpsn and ypsn can be grid-vectors, while upsn and vpsn are grid-arrays
%
% But these can also all be vectors, or all grid arrays. 
%
%%

if nargin > 2
    if isvector(xpsn) && size(xpsn,1) == size(upsn,1) && size(ypsn,1) == size(upsn,2)

        [xpsn,ypsn]=ndgrid(xpsn,ypsn);

    end
end


[lat,lon]=psn2ll(xpsn,ypsn);

[x,y]=LatLon2xyIceland(lat,lon); 

if nargin ==2
    u=[];
    v=[];
    return
end

%% rotating velocities...

% I'm guessing that Easting means in the direction of from (x0,y0) to (x0+d,y0) where d is some distance
%
% In the original polar  polar stereo graphic system, a vector point along the x axis will go from (xps0,yps0) to (xps0+dx,yps0) 
%
% If I find those same coordinates in the Icelandic coordinate system, the vector will, in general, not point along the
% x-axis of the Icelandic system, but be at some angle with respect to the Icelandic x axis. From this I can find the
% angle that I must rotate the vectors by. This is approximate, but should be very close. Possibly a more accurate
% approach would be to do this at each location of the velocity vectors anew.
%
%

[x0ps,y0ps]=ll2psn(65.0,-19.0) ;  % These are the psn coordinates of the lat/lon centre of the Icelandic coordinate system



d=1e3; % some small perturbation 
[lat0,lon0]=psn2ll(x0ps,y0ps);
[lat0,lon1]=psn2ll(x0ps+d,y0ps);

[x0Iceland,y0Iceland]=LatLon2xyIceland(lat0,lon0);
[x1Iceland,y1Iceland]=LatLon2xyIceland(lat0,lon1);

dd=([x1Iceland y1Iceland] - [x0Iceland y0Iceland]) ; 
dx=dd(1) ; dy=dd(2)  ;
theta=atan2d(dy,dx)  ;
theta=-theta;

u=upsn*cosd(theta)+vpsn*sind(theta);
v=-upsn*sind(theta)+vpsn*cosd(theta);


end