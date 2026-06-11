





function [lat,lon] = psn2ll(x,y,varargin)
% psn2ll transforms north polar stereographic map coordinates to geographic 
% lat,lon coordinates. This is an Antarctic Mapping Tools' ps2ll function, which 
% was derived from Andy Bliss' polarstereo_inv function. This function does NOT 
% require Matlab's Mapping Toolbox. 
% 
%% Syntax
% 
% [lat,lon] = psn2ll(x,y) 
% [lat,lon] = psn2ll(x,y,'TrueLat',ReferenceLatitude) 
% [lat,lon] = psn2ll(x,y,'EarthRadius',RadiusInMeters) 
% [lat,lon] = psn2ll(x,y,'Eccentricity',EarthsMisshapenness) 
% [lat,lon] = psn2ll(x,y,'meridian',MeridianInDegrees) 
% 
%% Description 
% 
% [lat,lon] = psn2ll(x,y) transforms polar stereographic x,y coordinates (re: 
% 70 N) to geographic lat/lon. Inputs x and y  can be scalar, vector, or
% matrices of equal size. 
% 
% [lat,lon] = psn2ll(x,y,'TrueLat',ReferenceLatitude) secifies a reference
% latitude of true scale in degrees; also known as the standard parallel.
% Default is 70 N. 
% 
% [lat,lon] = psn2ll(x,y,'EarthRadius',RadiusInMeters) specifies Earth's
% radius in meters. Default is 6378137.0 m, WGS84.
% 
% [lat,lon] = psn2ll(x,y,'Eccentricity',EarthsMisshapenness) specifies
% Earth's eccentricity or misshappenness.  Default values is 0.08181919. 
% 
% [lat,lon] = psn2ll(x,y,'meridian',MeridianInDegrees) specifies the meridian in 
% degrees along the positive Y axis of the map. Default value is -45.
%
%% Example: 
% Petermann Glacier is located at (80.75 N, 60.75 W). Get the north polar 
% stereographic coordinates in meters: 
% 
%    [xi,yi] = ll2psn(80.75,-60.75)
%    x =
%        -272559.27
%    y =
%        -966422.28
% 
% And now convert back to lat,lon coordinates: 
% 
%    [lat,lon] = psn2ll(x,y) 
%    lat =
%             80.75
%    lon =
%            -60.75
% 
%% Citing Antarctic Mapping Tools
% This function was originally developed for Antarctic Mapping Tools for Matlab. If it's useful for you,
% please cite our paper: 
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences. 104 (2017) pp.151-157. 
% http://dx.doi.org/10.1016/j.cageo.2016.08.003
% 
% @article{amt,
%   title={{Antarctic Mapping Tools for \textsc{Matlab}}},
%   author={Greene, Chad A and Gwyther, David E and Blankenship, Donald D},
%   journal={Computers \& Geosciences},
%   year={2017},
%   volume={104},
%   pages={151--157},
%   publisher={Elsevier}, 
%   doi={10.1016/j.cageo.2016.08.003}, 
%   url={http://www.sciencedirect.com/science/article/pii/S0098300416302163}
% }
%   
%% Futher Reading
%   
% Equations from: Map Projections - A Working manual - by J.P. Snyder. 1987 
% http://kartoweb.itc.nl/geometrics/Publications/Map%20Projections%20-%20A%20Working%20manual%20-%20by%20J.P.%20Snyder.pdf
% See the section on Polar Stereographic, with a south polar aspect and
% known phi_c not at the pole.
%
% WGS84 - radius: 6378137.0 eccentricity: 0.08181919
%   in Matlab: axes2ecc(6378137.0, 6356752.3142)
% Hughes ellipsoid - radius: 6378.273 km eccentricity: 0.081816153
%   Used for SSM/I  http://nsidc.org/data/polar_stereo/ps_grids.html
% International ellipsoid (following Snyder) - radius: 6378388.0 eccentricity: 0.0819919 
% 
% See also: ll2psn, PROJINV, PROJFWD, MINVTRAN, MFWDTRAN, and ROTATEM.
%% Input checks: 
assert(nargin>1,'The psn2ll function requires at least two inputs: mapx and mapy.')
%% Set defaults: 
phi_c = 70;   % standard parallel - this is different from Andy Bliss' function, which uses -70! 
a = 6378137.0; % radius of ellipsoid, WGS84
e = 0.08181919;% eccentricity, WGS84
lambda_0 = -45;  % meridian along positive Y axis
%% Parse user inputs: 
tmp = strncmpi(varargin,'true',4)|strncmpi(varargin,'lat',3)|strncmpi(varargin,'ref',3)|...
    strncmpi(varargin,'earthrad',8)|strncmpi(varargin,'rad',3)|...
    strncmpi(varargin,'ecc',3)|strncmpi(varargin,'PosYLon',4)|...
    strncmpi(varargin,'lon',3)|strncmpi(varargin,'merid',5); 
assert(sum(tmp)==(nargin-2)/2,'There seems to be at least one invalid input string. Are you trying to declare options that do not exist?')
tmp = strncmpi(varargin,'true',4)|strncmpi(varargin,'lat',3)|strncmpi(varargin,'ref',3); 
if any(tmp)
    phi_c = varargin{find(tmp)+1}; 
    assert(isscalar(phi_c)==1,'True lat must be a scalar.')
end
tmp = strncmpi(varargin,'earthrad',8)|strncmpi(varargin,'rad',3); 
if any(tmp)
    a = varargin{find(tmp)+1}; 
    assert(isscalar(a)==1,'Earth radius must be a scalar.')
    assert(a > 7e+3,'Earth radius should be something like 6378137 in meters.')
end
tmp = strncmpi(varargin,'ecc',3); 
if any(tmp)
    e = varargin{find(tmp)+1}; 
    assert(isscalar(e)==1,'Earth eccentricity must be a scalar.')
    assert(e>0 & e<1,'Earth eccentricity does not seem like a logical value.')
end
tmp = strncmpi(varargin,'PosYLon',4)|strncmpi(varargin,'lon',3)|...
    strncmpi(varargin,'merid',5); 
if any(tmp)
    lambda_0 = varargin{find(tmp)+1}; 
    assert(isscalar(lambda_0)==1,'PosYLon must be a scalar.')
    assert(lambda_0>=-180 & lambda_0<=360,'PsosYLon does not seem like a logical value.')
end
%% Transform: 
%convert to radians
phi_c=phi_c*pi/180;
lambda_0=lambda_0*pi/180;
%this is not commented very well. See Snyder for details.
t_c=tan(pi/4 - phi_c/2)/((1-e*sin(phi_c))/(1+e*sin(phi_c)))^(e/2);
m_c=cos(phi_c)/sqrt(1-e^2*(sin(phi_c))^2);
rho=sqrt(x.^2+y.^2); 
t=rho*t_c/(a*m_c);
%find phi with a series instead of iterating.
chi=pi/2 - 2 * atan(t);
phi=chi+(e^2/2 + 5*e^4/24 + e^6/12 + 13*e^8/360)*sin(2*chi)...
    + (7*e^4/48 + 29*e^6/240 + 811*e^8/11520)*sin(4*chi)...
    + (7*e^6/120+81*e^8/1120)*sin(6*chi)...
    + (4279*e^8/161280)*sin(8*chi);
lambda=lambda_0 + atan2(x,-y);
%correct the signs and phasing
lambda=mod(lambda+pi,2*pi)-pi; %want longitude in the range -pi to pi
%convert back to degrees
lat=phi*180/pi;
lon=lambda*180/pi;
%% Make two-column format if user requested fewer than two outputs: 
if nargout<2
    lat = [lat(:) lon(:)]; 
end
