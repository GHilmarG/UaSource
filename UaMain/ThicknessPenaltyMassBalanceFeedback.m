
function [aPenalty1,daPenaltydh1]=ThicknessPenaltyMassBalanceFeedback(CtrlVar,hint)

%%
% Calculates additional mass-balance term based on if ice thickness is below min ice thickness. This can be thought of as a
% penalty term. It is only applied if the ice thickness is below:
%
%   hmin=CtrlVar.ThickMin ;
% 
% This option is used in the -uvh- and the -h- solvers, and can be activated by setting:
%
%    CtrlVar.ThicknessPenalty=1;  
%
% The additional implicit mass-balance term has the form
%
% $$a^{\star} = a_1 (h-h_{\min}) + a_2 (h-h_{\min})^2 + a_3 (h-h_{\min})^3 $$
% 
% where
% 
% $$h < h_{\min}$$
% 
% and where:
%
%  a_1 = CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffLin a_2 = CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffQuad; a_3
%  = CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffCubic;
%
% Note: $a_1$ and $a_3$ need to be negative and $a_2$ positive, however, this is check internally so actually the sign on
% input is immaterial.
%
% This fictitious mass-balance term will be either positive or negative depending on whether the ice thickness is
% below or above that minimum ice thickness.
%
% For this term to be positive or negative depending on ice thickness with respect to the desired thickness, the functions
% are odd function of ice thickness and first and third power are allowed (but not second power).
%
%
%%

hmin=CtrlVar.ThickMin ;

% make sure the signs are correct.
a1= -abs(CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffLin);
a2= +abs(CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffQuad);
a3= -abs(CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffCubic);


PenaltyMask1=hint<hmin ;

% if thickness too small, then (hint-hmin) < 0, and ab > 0, provided a1 and a3 are negative

aPenalty1 =PenaltyMask1.* ( a1*(hint-hmin)+a2*(hint-hmin).^2 + a3*(hint-hmin).^3) ;
daPenaltydh1=PenaltyMask1.*(a1+2*a2*(hint-hmin) +3*a3*(hint-hmin).^2) ;


if CtrlVar.InfoLevelThickMin >= 1
    ThicknessLessThanZero=hint<0 ;
    if any(ThicknessLessThanZero)
        fprintf("\t Some hint negative at integration point %i with min(hint)=%g. \t Number of neg integration point thicknesses: %i \n",Iint,min(hint),numel(find(ThicknessLessThanZero)))
    end
end


end
