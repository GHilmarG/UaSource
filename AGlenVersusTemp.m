
function [AGlen,B]=AGlenVersusTemp(T)

    % Gives A as a function of temperature (degrees Celcius) in the units a^{-1} kPa^{-3}
    
    T=T+273.15;
    
    %      a0= 5.3d-15 % s-1 kPa-3
    % a0= 5.3d-24 % s-1 Pa-3 = m+3 s+5 kg-3
    a0 = 5.3e-15*365.2422*24.*60.*60. ; % a-1 kPa-3
    %       a0= 1.699d-61 ! m+3 a+5 kg-3
    
    % dimentionless temperature factor
    fa=1.2478766e-39 * exp(0.32769 * T) + 1.9463011e-10 * exp(0.07205 * T) ; %; %Smith & Morland 1982 ,  Keine Einheiten; fa(273.15)=1.0
    
    AGlen=a0 * fa  ;  % 
    

    rhog=917*9.81/1000 ;   % kPa/m
    B=0.5*AGlen*(rhog)^3 ;  %  m^{-3} a^{-1} 


%  eta=10^15 Pa s = 10^13 kPa s = 10^13 /(365.25*24*60*60) kPa yr


end
