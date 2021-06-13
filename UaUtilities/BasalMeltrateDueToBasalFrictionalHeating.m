function ab=BasalMeltrateDueToBasalFrictionalHeating(ub,taub,rho,L)

%%

%  [ub] = m/yr
%  [taub] = kPa = kJ/m^3
%
%  [u tau ]=  m/s * Pa = m / s * J/m^3 = J s /m^2
%
% Enthalpy of fusion = L = [J/gramm = kJ/kg]
% [rho] = kg/m^3
%
% [rho L]=  kg/m^3 kJ/kg = kJ/m^3
%
%  u tau/(L *rho) = m/yr  kP  kg/kJ m^3/kg
%                  =m/yr  kJ/m^3 m^3/kJ
%                  =m/yr
%
% Calculating basal melt rates due to frictional heating in Ua is simple: 
%
%   tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,F.ub,F.vb,F.C,F.m,F.GF);
%   L=333.44 ; % Enthalpy of fusion = L = [J/gramm = kJ/kg]
%   ab=(tbx.*F.ub+tby.*F.vb)./(F.rho.*L);
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
