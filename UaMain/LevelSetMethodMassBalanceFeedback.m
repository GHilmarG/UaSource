
function [abLSF,dabLSFdh]=LevelSetMethodMassBalanceFeedback(CtrlVar,LM,hint)

%%
% Calculates additional mass-balance forcing to cause ice thickness downstream of the calving front to be close to:
%
%   CtrlVar.LevelSetMinIceThickness;
%
% LM is the level-set mask which must be 1 downstream of calving fronts, and zero elsewhere. 
%
%%

narginchk(3,3)

a1= -abs(CtrlVar.LevelSetMethodMassBalanceFeedbackCoeffLin);
a3= -abs(CtrlVar.LevelSetMethodMassBalanceFeedbackCoeffCubic);



hmin=CtrlVar.LevelSetMinIceThickness;

abLSF =LM.* ( a1*(hint-hmin)+a3*(hint-hmin).^3) ;
dabLSFdh=LM.*(a1+3*a3*(hint-hmin).^2) ;


% Note: This is a change done in Sept 2026 where I now divide by dt.
abLSF=abLSF/(CtrlVar.dt+eps(CtrlVar.dt)) ;
dabLSFdh=dabLSFdh/(CtrlVar.dt+eps(CtrlVar.dt)) ;



end