
function [aPenalty1,daPenaltydh1]=ThicknessPenaltyMassBalanceFeedback(CtrlVar,hint)

persistent WarnedAboutPenaltyWidth

%%
% Calculates additional mass-balance term based on if ice thickness is below min ice thickness. This can be thought of as a
% penalty term. It is only applied if the ice thickness is below hmin, where for the polynomial and the exponential options
%
%   hmin=2*CtrlVar.ThickMin ;
%
% and for the SoftPlus option (the default option, see below)
%
%   hmin=CtrlVar.ThickMin+delta,  delta=max(CtrlVar.ThickMin,CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.deltaAbs)
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
% Setting
%
%   CtrlVar.ThicknessPenaltyMassBalanceFeedbackFunction="softplus";
%
% adds a smooth linear trend
%
% $$ a=K \; f(k,-h,-h_{\min}) $$
%
%
% where
%
% $$ f(k,x,x_0) =\mathrm{SoftPlus}(x) = \frac{1}{2k} \, \ln \left ( 1+e^{2k(x-x_0)} \right ) $$
%
% The function $f$ is the SoftPlus function. It is a smooth version of a function that is zero for $x < x_0$ with a smoothness determined
% by  $k$ where $h$ has the units of inverse $x$.
%
% So, for example,
%
% $$a=K \;  f(1/(2 l),-h,-h_{\min}) $$
%
% gives $a$  that is zero for $h>h_{\min}$ and approximately  equal to $-K (h-h_{\min})$ for $h<h_{\min}$. This gives a positive mass
% balance for $h$ smaller than $h_{\min}$ and zero if $h$ is greater than $h_{\min}$, with a smoothing distance scale $l$
%
%   K= CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.K;
%   l= CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.l;
%
% K is a rate (units 1/time). If l is not set (default l=NaN), it is calculated at runtime as l=lRel*delta. The default values are:
%
%  CtrlVar.ThicknessPenaltyMassBalanceFeedbackFunction="softplus";
%  CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.K=100;     % 1/time
%  CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.l=NaN;     % ie l=lRel*delta
%  CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.lRel=0.1;
%  CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.deltaAbs=0.1;
%
% At h=ThickMin the penalty is then in its linear range and supplies a mass-balance rate of about K*delta. A node is held above
% ThickMin by the penalty alone if theta*K*delta > |a|max (theta=CtrlVar.theta), otherwise the active set takes over.
%
%
%
%
%%

switch lower(CtrlVar.ThicknessPenaltyMassBalanceFeedbackFunction)

    case "polynomial"

        %% Polynomial barrier
        hmin=2*CtrlVar.ThickMin ;

        % make sure the signs are correct.
        a1= -abs(CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffLin);
        a2= +abs(CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffQuad);
        a3= -abs(CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffCubic);


        %PenaltyMask1=hint<hmin ;
        k=10000/hmin;
        PenaltyMask1 = HeavisideApprox(k,hmin,hint) ;
        dPenaltyMask1dh = -DiracDelta(k,hmin,hint) ;

        % if thickness too small, then (hint-hmin) < 0, and ab > 0, provided a1 and a3 are negative

        aPenalty1 = PenaltyMask1.* ( a1*(hint-hmin)+a2*(hint-hmin).^2 + a3*(hint-hmin).^3) ;
        daPenaltydh1=PenaltyMask1.*(a1+2*a2*(hint-hmin) +3*a3*(hint-hmin).^2) +  dPenaltyMask1dh .* ( a1*(hint-hmin)+a2*(hint-hmin).^2 + a3*(hint-hmin).^3) ;


    case "exponential"

        %% exponential barrier
        K= CtrlVar.ThicknessPenaltyMassBalanceFeedbackExponential.K;
        l= CtrlVar.ThicknessPenaltyMassBalanceFeedbackExponential.l;
        hmin=2*CtrlVar.ThickMin ;
        aPenalty1=K*exp(-(hint-hmin)/l);
        daPenaltydh1=-K*exp(-(hint-hmin)/l)/l;


    case "softplus"
        %% Softplus

        % (9 Oct 2026) The penalty is a mass-balance rate, ie K has units of 1/time. Between 31 Aug and 9 Oct 2026, K was scaled with
        % 1/dt, making the penalty a per-time-step term that changes the thickness by a dt-independent amount in each time step. As a
        % consequence the discrete solution had no limit for dt->0, the number of Newton iterations did not decrease for small dt,
        % and the automated time stepping could drive dt to very small values. This scaling has been removed. When dt is small, the
        % active set (thickness constraints) does the work that the per-time-step penalty was intended to do. A node is held above
        % ThickMin by the penalty alone if theta*K*delta > |a|max, see Ua2D_DefaultParameters.

        % Modification (8 Oct 2026): penalty centred a distance delta above ThickMin,
        % delta=max(ThickMin,deltaAbs). For ThickMin>=deltaAbs this is the previous hmin=2*ThickMin,
        % but it does not collapse to hmin=0 for ThickMin=0.
        if isfield(CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus,"deltaAbs") ; deltaAbs=CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.deltaAbs ; else ; deltaAbs=0.1 ; end
        hmin=CtrlVar.ThickMin+max(CtrlVar.ThickMin,deltaAbs) ;

        K= CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.K ;   % a rate, units 1/time
        % (9 Oct 2026) The smoothing distance l is calculated at runtime from delta, l=lRel*delta, unless it has been set explicitly
        % (default l=NaN). With lRel=0.1 the penalty is in its linear range at h=ThickMin, ie (hmin-ThickMin)/l=10, and it only affects
        % the thickness in a narrow band above ThickMin.
        delta=max(CtrlVar.ThickMin,deltaAbs) ;
        l= CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.l;
        if isempty(l) || isnan(l)
            lRel=0.1 ;
            if isfield(CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus,"lRel") ; lRel=CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.lRel ; end
            l=lRel*delta ;
        elseif l>delta/2 && isempty(WarnedAboutPenaltyWidth)
            warning("ThicknessPenaltyMassBalanceFeedback:Width",...
                "CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.l=%g is larger than delta/2=%g, where delta=max(ThickMin,deltaAbs). The penalty is then not in its linear range at h=ThickMin, and its capacity differs from K*delta. Consider l=0.1*delta, or leave l unset (NaN) so that it is calculated as lRel*delta.",l,delta/2)
            WarnedAboutPenaltyWidth=true;
        end

        k=1/(2*l);

        [aPlus,daPlusdh] = SoftPlus(k,-hint,-hmin);
        aPenalty1=K*aPlus;
        daPenaltydh1=-K*daPlusdh ; % don't forget the outer derivative, because the input to SoftPlus is -hint and not +hint


    otherwise

        error("ThicknessPenaltyMassBalanceFeedback:CaseNotFound","Case not found")

end


if CtrlVar.InfoLevelThickMin >= 10
    ThicknessLessThanZero=hint<0 ;
    if any(ThicknessLessThanZero)
        fprintf("\t Some hint negative at integration point with min(hint)=%g. \t Number of neg integration point thicknesses: %i \n",min(hint),numel(find(ThicknessLessThanZero)))
    end

    if CtrlVar.InfoLevelThickMin >= 10
        %%
        Fig=FindOrCreateFigure("a penalty versus ice thickness") ; clf(Fig)

        yyaxis left ;
        plot(hint,aPenalty1,".b",DisplayName="$a_{\mathrm{penalty}}$")
        ylabel("a Penalty")
        hold on ;

        yyaxis right ;
        plot(hint,daPenaltydh1,".r",DisplayName="$da_{\mathrm{penalty}}/dh$")
        ylabel("da/dh")

        xrange=[min(0,min(hint)) 5*hmin];
        if xrange(1)==xrange(2)
            xrange(1)=hmin-1;
            xrange(2)=hmin+1; 
        end
        xlim(xrange) ;
        xlabel("hint") ;
        title("a penalty") ;
        xline(hmin,"--k","hmin",DisplayName="$h_{\mathrm{min}$")
        %ylim([-K/l 0])


        switch lower(CtrlVar.ThicknessPenaltyMassBalanceFeedbackFunction)

            case "softplus"
                hExample=linspace(xrange(1),xrange(2)) ;

                % k=1/l ;
                [aPlusExample,daPlusdhExample] = SoftPlus(k,-hExample,-hmin);
                aPlusExample=K*aPlusExample;
                daPlusdhExample=-K*daPlusdhExample ; %
        end
        hold on ; yyaxis left ; plot(hExample,aPlusExample,"--",DisplayName="$a_{\mathrm{penalty}}$ curve")
        hold on ; yyaxis right ; plot(hExample,daPlusdhExample,"--",DisplayName="$da_{\mathrm{penalty}}/dh$ curve")
        lg=legend(Interpreter="latex"); 
        %%
    end

end

% FindOrCreateFigure("penalty mask versus ice thickness") ; plot(hint,PenaltyMask1,".k") ; xlim([hmin-5/k hmin+5/k]) ; title("Penalty mask")
% FindOrCreateFigure("a penalty versus ice thickness") ; plot(hint,aPenalty1,".k") ; xlim([0 2*hmin]) ; title("a penalty") ; xline(hmin,"--k")

% numel(find(hint<hmin))



end
