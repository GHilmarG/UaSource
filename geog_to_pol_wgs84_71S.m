%Transforms from geodetic coordinates (lat,lon) to 
%Southern Hemisphere polar stereographic coordinates (x,y)
%with latitude of true scale at 71S.
%Reference:
%Snyder, John P. Map Projections Used by the United States Geological Survey-2nd edition.
%U.S.G.S. Bulletin No. 1532. Washington D.C.: U.S. Government Printing Office, 1983. 
function [x,y]=geog_to_pol_wgs84_71S(lat,lon,phi_c,lon0)

%conversion factor from degrees to radians
dtorad=pi/180;


if nargin < 3
%latitude of true scale
    phi_c=-71;
end

if nargin<4
    lon0=0;
end

%WGS84 ellipsoid
a=6378137;
f=1/298.257223563;
e2=2*f-f^2;
e=sqrt(e2);

%Snyder equations 12-15 and 13-9 evaluated at latitude of true scale
m_c=cos(phi_c*dtorad)/sqrt(1-e2*(sin(phi_c*dtorad))^2);  
t_c=tan(pi/4+phi_c*dtorad/2)/ ( (1+e*sin(phi_c*dtorad))/(1-e*sin(phi_c*dtorad)) )^(e/2);

%Snyder equation 13-9 evaluated at specified latitudes;
t=tan(pi/4+lat*dtorad/2)./ ( (1+e*sin(lat*dtorad))./(1-e*sin(lat*dtorad)) ).^(e/2);

%Snyder equation 17-34
rho=a*m_c*t/t_c;

%Snyder equations 17-30 and 17-31
x=-rho.*sin((lon0-lon)*dtorad);
y=rho.*cos((lon0-lon)*dtorad);


