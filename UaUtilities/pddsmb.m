function [smb] = pddsmb(h,T,varargin)
%
% calculates smb as a function of temperature using a simple pdd style
% model. 
%
% smb = pddsmb(z,Tsealevel [,param,value]) 
%
% Parameters:
%  * Tlapse: Temperature lapse rate [-7e-3]
%  * Tsigma: Temperature standard deviation [7.5] (should include both and hi and hifreq noise)
%  * PDDfactor: [8e-3];
%  * P0: precipitation at sea level. Will be reduced by 9% per degree
%  
%
% Aslak Grinsted

p=inputParser;
p.FunctionName='pddsmb';
try %older versions of matlab did not support partialmatching.
    p.CaseSensitive=false;
    p.StructExpand=true;
    p.PartialMatching=true;
catch
end
p.addOptional('Tlapse',-7e-3,@isnumeric); %-7e-3 temperature lapse rate (http://www.igsoc.org:8080/journal/55/189/t08J033.pdf)
p.addParameter('Tsigma',7.5,@isnumeric); %this is a sigma that should capture both annual cycle and hi freq noise
                                         %standard deviation of temperature judged from: http://www.arctic.noaa.gov/report12/images-terrcryo/g-fig5.16.jpg
p.addParameter('PDDfactor',8e-3,@isnumeric); %8 is pretty good according to: http://www.igsoc.org:8080/journal/41/137/igs_journal_vol41_issue137_pg153-160.pdf
p.addParameter('P0',2.5,@isnumeric); %precip at T=0
p.parse(varargin{:});
R=p.Results;


% ------------- surface temperature
Tsurf=T+R.Tlapse*h;

%------------- ABLATION --------------

fPDD=normalpdd(Tsurf,R.Tsigma); %standard deviation of temperature judged from: http://www.arctic.noaa.gov/report12/images-terrcryo/g-fig5.16.jpg
abl=R.PDDfactor*fPDD*365.25;

%------------ ACCUMULATION -------------

P=R.P0*exp(0.09*(Tsurf)); %the 9% is based on fig a2 here http://www.clim-past.net/9/1029/2013/cp-9-1029-2013.pdf

%snow fraction
Pcold=normcdf(1,Tsurf,R.Tsigma);
acc=P.*Pcold;
smb=acc-abl;







function fPDD=normalpdd(Tmu,sigma)
% PDD fraction in a normally distributed climate
%
% Calculates the probability weighted positive degree days for normally distributed temperatures:
%  fPDD=int(normpdf(T,Tmu,sigma)*T,'T',0,inf)
%
% Usage: 
%  fPDD=normalpdd(Tmu,sigma)
%
% Example:
%   T_june = 3; %degC
%   sigma_june = 7; %degC  =std(daily_average_temps_june)
%   PDD_june=30*normalpdd(T_june,sigma_june)  % units: degC*days
%   PDD_factor = 8e-3 %m/day/degC   <- different factors exist in literature. 
%   melt_june = PDD_factor*PDD_june % meters of melt.
%
% reference: Braithwaite 1985. 
%
% Aslak Grinsted 2012


zz=Tmu./(sqrt(2)*sigma);
fPDD=(sigma/(sqrt(2*pi)))./exp(zz.^2) + Tmu.*(0.5+0.5*erf(zz));


