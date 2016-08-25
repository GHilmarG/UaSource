
function [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)


%%
%  Define mass-balance distributions along upper and lower ice surfaces:
%
%    [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
%    [UserVar,as,ab,dasdh,dabdh]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
%  
%   as     :  upper surface mass balance ab     :  lower surface mass balance
%   dasdh  :  gradient of upper surface mass balance with respect to ice
%   thickness dabdh  :  gradient of lower surface mass balance with respect to
%   ice thickness
%
% The mass balance is given in the units  in units distance/time, for example
% m/yr:
% 
%
% If mass-balance geometry feedback is included in the solver, define mass
% balance as:
%    [UserVar,as,ab,dasdh,dabdh]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
% 
% If this feedback is not included in the solver there is no need to return
% either dasdh nor dabdh:
%
% Note that the difference between including and not including the mass-balance
% feedback in the solver is a tecknical one. If the feeback is included in the
% solver, then the Newton-Raphson system is modified accordingly and these
% derivatives included in the left-hand side. This speeds up the solution of the
% non-linear system and increases the accuracy. At the same time, not specifying
% this feedback will not necessarly make the solution incorrect.   Even if the
% mass-balance feedback is not treated implicity, the mass balance is always
% calculated on the basis of the current geometry and this feedback therefore
% (explicitly) included.
%
%%

as=s*0;
ab=s*0;

end

