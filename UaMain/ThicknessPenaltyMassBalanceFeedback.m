
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
%   a_1 = CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffLin 
%   a_2 = CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffQuad; 
%   a_3= CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffCubic;
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

switch lower(CtrlVar.ThicknessPenaltyMassBalanceFeedbackFunction)

    case "polynomical"

        %% Polynomial barrier
        hmin=2*CtrlVar.ThickMin ;

        % make sure the signs are correct.
        a1= -abs(CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffLin);
        a2= +abs(CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffQuad);
        a3= -abs(CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffCubic);


        %PenaltyMask1=hint<hmin ;
        k=10000/hmin;
        PenaltyMask1 = HeavisideApprox(k,hmin,hint) ;
        dPenaltyMask1dh = DiracDelta(k,hmin,hint) ;

        % if thickness too small, then (hint-hmin) < 0, and ab > 0, provided a1 and a3 are negative

        aPenalty1 = PenaltyMask1.* ( a1*(hint-hmin)+a2*(hint-hmin).^2 + a3*(hint-hmin).^3) ;
        daPenaltydh1=PenaltyMask1.*(a1+2*a2*(hint-hmin) +3*a3*(hint-hmin).^2) +  dPenaltyMask1dh .* ( a1*(hint-hmin)+a2*(hint-hmin).^2 + a3*(hint-hmin).^3) ;


    case "exponential"

        %% exponential barrier
        K= CtrlVar.ThicknessPenaltyMassBalanceFeedbackExponential.K;
        l= CtrlVar.ThicknessPenaltyMassBalanceFeedbackExponential.l;
        hmin=CtrlVar.ThickMin ; 
        %K=10; l=hmin/10 ;
        aPenalty1=K*exp(-(hint-hmin)/l);
        daPenaltydh1=-K*exp(-(hint-hmin))/l;


    case "softplus"
        %% Softplus

        K= CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.K;
        l= CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.l;

        hmin=CtrlVar.ThickMin ; 
        %K=2000 ; l=hmin/10;
        E=exp(-(hint-hmin)/l);
        SoftPlus=K*log(1+E);
        dSoftPlusdh=(-K/l)./ (1./E+1);
        aPenalty1=SoftPlus;
        daPenaltydh1=dSoftPlusdh;

    otherwise

        error("ThicknessPenaltyMassBalanceFeedback:CaseNotFound","Case not found")

end


if CtrlVar.InfoLevelThickMin >= 1
    ThicknessLessThanZero=hint<0 ;
    if any(ThicknessLessThanZero)
        fprintf("\t Some hint negative at integration point %i with min(hint)=%g. \t Number of neg integration point thicknesses: %i \n",min(hint),numel(find(ThicknessLessThanZero)))
    end

    if CtrlVar.InfoLevelThickMin >= 10
        Fig=FindOrCreateFigure("a penalty versus ice thickness") ; clf(Fig)

        yyaxis left ;
        plot(hint,aPenalty1,".b") ;
        ylabel("a Penalty")
        hold on ;

        yyaxis right ;
        plot(hint,daPenaltydh1,".r") ; ylabel("da/dh")

        xlim([0 5*hmin]) ;
        xlabel("hint") ;
        title("a penalty") ;
        xline(hmin,"--k","hmin") ;
        ylim([-K/l 0])

    end

end

% FindOrCreateFigure("penalty mask versus ice thickness") ; plot(hint,PenaltyMask1,".k") ; xlim([hmin-5/k hmin+5/k]) ; title("Penalty mask")
% FindOrCreateFigure("a penalty versus ice thickness") ; plot(hint,aPenalty1,".k") ; xlim([0 2*hmin]) ; title("a penalty") ; xline(hmin,"--k")

% numel(find(hint<hmin))



end
