
function [AGlen,B]=AGlenVersusTemp(T)


%%
% Gives A as a function of temperature (degrees Celsius) in the units a^{-1} kPa^{-3}
%
%   $$\dot{\epsilon}=A \tau^m$$
%
%
% Example:
%
% AGlen=AGlenVersusTemp(0)   ; % Gives A in Glen's flow law for n=3 and T=0 degrees Celsius in the units kPa^{-3} yr^{-1}
%
% Based on Smith & Morland, 1982
%
%%


T=T+273.15;

%       a0= 5.3d-15 ; % s-1 kPa-3
%       a0= 5.3d-24 ; % s-1 Pa-3 = m+3 s+5 kg-3
a0 = 5.3e-15*365.2422*24.*60.*60. ; % a-1 kPa-3
%       a0= 1.699d-61 ; % m+3 a+5 kg-3

% dimensionless temperature factor
fa=1.2478766e-39 * exp(0.32769 * T) + 1.9463011e-10 * exp(0.07205 * T) ; %; %Smith & Morland 1982 ,  no units, scaled such that:  fa(273.15)=1.0

AGlen=a0 * fa  ;  %

if nargout>1
    rhog=917*9.81/1000 ;   % kPa/m
    B=0.5*AGlen*(rhog)^3 ;  %  m^{-3} a^{-1}, some publications use this definition 
else
    B=[];
end


end
