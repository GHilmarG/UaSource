function ab=BasalMeltrateDueToBasalFrictionalHeating(ub,taub,rho,L)

%%

%  [ub] = m/yr
%  [taub] = kPa
%
%  [u tau ]=  m/s * Pa = m / s * J/m^3 = J s /m^2
%
% Enthalpy of fusion = L = [J/gramm = kJ/kg]
% [rho] = kg/m^3
%
%  u tau/(L *rho) = m/yr  kP  kg/kJ m^3/kg
%                  =m/yr  kJ/m^3 m^3/kJ
%                  =m/yr
%

%
%%

if nargin <3
    rho=910 ;
    L=333.44;
elseif nargin <4
    L=333.44;
end


ab=ub*taub/(L*rho); 




end
