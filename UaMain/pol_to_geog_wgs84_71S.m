function [lat,lon]=pol_to_geog_wgs84_71S(x,y)

%Transforms from Southern Hemisphere polar
%stereographic coordinates (x,y)
%with latitude of true scale at 71S
%to geodetic coordinates (lat,lon).
%Reference:
%Snyder, John P. Map Projections Used by the United States Geological Survey-2nd edition.
%U.S.G.S. Bulletin No. 1532. Washington D.C.: U.S. Government Printing Office, 1983.



%conversion factor from degrees to radians
dtorad=pi/180;
radtod=180/pi;

%latitude of true scale
phi_c=-71;
lon0=0;

%WGS84 ellipsoid
a=6378137;
f=1/298.257223563;
e2=2*f-f^2;
e=sqrt(e2);

%Snyder equations 12-15 and 13-9 evaluated at latitude of true scale
m_c=cos(phi_c*dtorad)/sqrt(1-e2*(sin(phi_c*dtorad))^2);
t_c=tan(pi/4+phi_c*dtorad/2)/ ( (1+e*sin(phi_c*dtorad))/(1-e*sin(phi_c*dtorad)) )^(e/2);

%Snyder equation 16-18
rho=sqrt(x.^2+y.^2);
t=rho*t_c/(a*m_c);

tol=1.0d-15;
phi_diff=tol;
%phi_old=-(pi/2-2*atan(t));
phi_old=-(pi/2-2*atan2(rho*t_c,a*m_c));

while (phi_diff >= tol);
    %phi_new = -(pi/2)+2*atan(t .* ( (1+e*sin(phi_old)) ./ (1-e*sin(phi_old)) ).^(e/2) );
    phi_new = -(pi/2)+2*atan2 (t .* (1+e*sin(phi_old)).^(e/2) , (1-e*sin(phi_old)).^(e/2) );
    phi_diff = max (max (abs(phi_new - phi_old) ) );
    phi_old=phi_new;
end

lat=phi_new*radtod;

lon=lon0-radtod*atan2(-x,y);
polrad=1;
pole=rho <= polrad;
lon(pole)=0;

% q2=find((y < 0) & (x >= 0));
% lon(q2)=lon(q2)+180;
% q3=find((y < 0) & (x < 0));
% lon(q3)=lon(q3)-180;

