

function [gammaMin,Jmin,info]=LineSearchWolfe(slope0,gammaStart,J0,J1,Func,info)

%%
%
%   [gammaMin,Jmin,info]=LineSearchWolfe(slope0,gammaStart,J0,J1,Func,info)
%
% Line search returning a step gammaMin that satisfies the strong Wolfe
% conditions. Intended for use with non-linear conjugated gradients, where the
% curvature condition must be tight for the conjugacy of successive search
% directions to be preserved.
%
% Writing  phi(gamma) = J(p0 + gamma*d)  the strong Wolfe conditions are
%
%   1) sufficient decrease (Armijo)   phi(gamma) <= phi(0) + c1 gamma phi'(0)
%   2) curvature                      |phi'(gamma)| <= c2 |phi'(0)|
%
% This is the bracketing/zoom algorithm of Nocedal & Wright, Algorithms 3.5 and
% 3.6, with safeguarded cubic and quadratic interpolation. It is the same
% approach as MINPACK dcsrch / More-Thuente.
%
%
%% Inputs
%
%   slope0        phi'(0) , the directional derivative at gamma=0. A SCALAR,
%                 and it must be negative, i.e. d must be a descent direction.
%   gammaStart    initial trial step
%   J0            phi(0)          , already known to the caller
%   J1            phi(gammaStart) , already known to the caller
%   Func          function handle, [J,slope]=Func(gamma) , where J=phi(gamma)
%                 and slope=phi'(gamma) . BOTH must be scalars.
%   info          (optional) structure of options, see below
%
%
%% Building the directional derivative
%
% slope is dJ/dgamma along the search direction d, NOT the parameter gradient
% dJ/dp . If the cost function evaluator returns the gradient with respect to p,
% as is usual, it has to be contracted with d before being passed in here:
%
%       phi'(gamma) = DJ(p0+gamma*d)[d]
%
% Which contraction is correct depends on which gradient is at hand:
%
%   raw l2 gradient g, as returned by the adjoint      phi' = g'*d
%   Riesz-mapped gradient dJdp = G\g                   phi' = dJdp'*(G*d)
%
% The two agree, since (G*dJdp)'*d = g'*d . Mixing them up rescales the slope
% and silently breaks both Wolfe conditions.
%
% At gamma=0 there is a shortcut. With s1=-dJdp the steepest-descent direction
% and d1 the search direction returned by NewConjugatedGradMetric,
%
%       slope0 = g'*d1 = -<s1,d1>_G
%
% which is minus the sufficient-descent quantity already formed inside that
% routine, so it can be taken straight from there rather than recomputed.
%
% Outputs
%
%   gammaMin      accepted step
%   Jmin          phi(gammaMin)
%   info          structure, with
%
%                   info.nFuncEvaluations   number of calls to Func made here.
%                                           The two evaluations giving J0 and J1
%                                           are NOT counted, having been made by
%                                           the caller.
%                   info.nGradEval          number of those calls for which the
%                                           slope was requested, i.e. the number
%                                           of gradient evaluations. See the note
%                                           on counting below.
%                   info.slopeMin           phi'(gammaMin)
%                   info.Armijo             true if condition 1) holds
%                   info.Curvature          true if condition 2) holds
%                   info.iExit              see below
%                   info.Message            text describing iExit
%
%
%% Options, all fields of info, set to defaults if absent or empty
%
%   info.c1                   1e-4    Armijo parameter
%   info.c2                   0.1     curvature parameter. 0.1 is the usual
%                                     choice for CG (Nocedal & Wright section
%                                     5.2); 0.9 is for quasi-Newton and is far
%                                     too loose to keep CG directions conjugate.
%                                     Must satisfy 0 < c1 < c2 < 1.
%   info.MaxFuncEvaluations   15
%   info.gammaMax             1e6*gammaStart   upper bound on the step. Tying
%                                     this to gammaStart is only a fallback: a
%                                     poor initial guess then becomes a hard
%                                     ceiling and the search exits with iExit=4
%                                     having satisfied Armijo but not the
%                                     curvature condition. Set it from the
%                                     problem where a genuine bound exists, e.g.
%                                     a trust-region radius or parameter bounds.
%   info.xtol                 1e-10   bracket collapse tolerance, relative
%   info.slopeStart           []      phi'(gammaStart) if the caller already has
%                                     it. Supplying it saves one evaluation.
%   info.InfoLevel            0
%
%
%% Exit flags
%
%    1   strong Wolfe conditions satisfied
%    2   evaluation limit reached, returning best point found, Armijo satisfied
%    3   bracket collapsed before the curvature condition was met. Usually means
%        phi is very flat, or noisy at the level of the requested accuracy.
%    4   gammaMax reached
%   -1   slope0 >= 0 , the direction is not a descent direction
%   -2   no point better than phi(0) was found
%
%
%% A warning on non-scalar inputs
%
% The inputs are validated because getting this wrong fails quietly rather than
% loudly. In MATLAB an if statement on an array is true only when every element
% is nonzero, so passing a gradient vector where a scalar slope is expected does
% not make the curvature test fire early, it makes it almost never fire, and the
% search then runs to MaxFuncEvaluations and exits with iExit=2. Worse, the
% descent check ~(slope0<0) would then pass unless every component of slope0
% were non-negative, so a genuine ascent direction could go undetected.
%
%
%% Notes
%
% The curvature condition cannot be tested without phi'(gamma), so Func must
% return the slope as its second output. In a PDE-constrained inverse problem
% the slope at a trial point costs one adjoint solve on top of the forward solve
% already performed to get J, so a [J,slope] evaluation is considerably cheaper
% than two J evaluations.
%
%
%% A note on counting: why nGradEval equals nFuncEvaluations
%
% Both counters are returned, but they are necessarily equal here, and it is
% worth being clear why rather than leaving it looking like an oversight.
%
% Whether the slope at a trial point is needed depends on which branch of the
% algorithm that point falls into, and that is only known after J has been
% computed there. In the zoom phase a point failing the Armijo or decrease test
% becomes the new bracket end, where the slope is used only to allow cubic
% rather than quadratic interpolation; a point passing it becomes the new low
% end, where the slope is needed both to test the curvature condition and to
% decide which way the bracket collapses. Since the branch cannot be predicted
% in advance, deferring the gradient would mean either calling Func twice at the
% same gamma, so paying for a second forward solve to save one adjoint solve, or
% abandoning cubic interpolation altogether. Neither is worth it, so the slope
% is always requested and the two counters coincide.
%
% They are reported separately all the same, so that a caller accumulating
% totals across a run does not have to know this, and so that the fields stay
% meaningful if the evaluation strategy is ever changed.
%
% Note also that the Armijo test at gammaStart uses the caller-supplied J1 and
% needs no gradient. When that test fails the routine drops straight into the
% zoom phase having called Func not at all, so the gradient at gammaStart is
% never computed in that case.
%
% Unlike backtracking, this routine may increase the step beyond gammaStart. A
% pure backtracking search can only ever satisfy the Armijo condition and gives
% no control at all over the curvature condition, which is why it is a poor
% match for CG.
%
%%

if nargin<6 || ~isstruct(info) ; info=struct ; end

info=SetDefault(info,'c1',1e-4) ;
info=SetDefault(info,'c2',0.1) ;
info=SetDefault(info,'MaxFuncEvaluations',15) ;
info=SetDefault(info,'gammaMax',1e6*gammaStart) ;
info=SetDefault(info,'xtol',1e-10) ;
info=SetDefault(info,'slopeStart',[]) ;
info=SetDefault(info,'InfoLevel',0) ;

c1=info.c1 ; c2=info.c2 ; gammaMax=info.gammaMax ;

if ~(0<c1 && c1<c2 && c2<1)
    error('LineSearchWolfe:BadParameters','Require 0<c1<c2<1, but c1=%g and c2=%g.',c1,c2)
end

% Validate that the scalars really are scalars. See the warning in the header:
% non-scalar slopes do not error of their own accord, they quietly change the
% meaning of every if statement below.
CheckScalar(slope0,'slope0') ;
CheckScalar(J0,'J0') ;
CheckScalar(J1,'J1') ;
CheckScalar(gammaStart,'gammaStart') ;
if ~isempty(info.slopeStart) ; CheckScalar(info.slopeStart,'info.slopeStart') ; end

if ~(gammaStart>0)
    error('LineSearchWolfe:BadGammaStart','gammaStart must be positive, but is %g.',gammaStart)
end

info.nFuncEvaluations=0 ;
info.nGradEval=0 ;
info.iExit=0 ;
info.Message='' ;

% best point seen so far, used if the search has to give up
gammaBest=0 ; JBest=J0 ; slopeBest=slope0 ;

if ~(slope0<0)
    gammaMin=0 ; Jmin=J0 ;
    info.iExit=-1 ;
    info.Message='slope0>=0, not a descent direction.' ;
    info.slopeMin=slope0 ; info.Armijo=false ; info.Curvature=false ;
    return
end

%% bracketing phase, Nocedal & Wright Algorithm 3.5

gammaPrev=0        ; JPrev=J0 ; slopePrev=slope0 ;
gammaCur=gammaStart ; JCur=J1  ; slopeCur=info.slopeStart ;

iBracket=0 ;

while true

    iBracket=iBracket+1 ;

    UpdateBest(gammaCur,JCur,slopeCur) ;

    % Armijo violated, or no decrease relative to the previous trial point:
    % a minimiser is bracketed by [gammaPrev,gammaCur]
    if JCur > J0+c1*gammaCur*slope0 || (iBracket>1 && JCur >= JPrev)
        [gammaMin,Jmin]=Zoom(gammaPrev,JPrev,slopePrev,gammaCur,JCur) ;
        break
    end

    if isempty(slopeCur)
        [JCur,slopeCur]=Evaluate(gammaCur) ;
        UpdateBest(gammaCur,JCur,slopeCur) ;
    end

    % curvature condition met, accept
    if abs(slopeCur) <= -c2*slope0
        gammaMin=gammaCur ; Jmin=JCur ;
        info.slopeMin=slopeCur ; info.iExit=1 ;
        info.Message='strong Wolfe conditions satisfied.' ;
        break
    end

    % slope has turned positive, minimiser is bracketed by [gammaCur,gammaPrev]
    if slopeCur >= 0
        [gammaMin,Jmin]=Zoom(gammaCur,JCur,slopeCur,gammaPrev,JPrev) ;
        break
    end

    if gammaCur >= gammaMax
        gammaMin=gammaCur ; Jmin=JCur ;
        info.slopeMin=slopeCur ; info.iExit=4 ;
        info.Message='gammaMax reached.' ;
        break
    end

    if info.nFuncEvaluations >= info.MaxFuncEvaluations
        [gammaMin,Jmin]=GiveUp ;
        break
    end

    % Still descending. Extrapolate using a cubic where that is sensible, but
    % never by less than the factor xtrapf, which is the MINPACK dcsrch value.
    % Without this floor a badly scaled gammaStart takes many doublings to
    % recover.
    xtrapf=4 ;
    gammaNew=CubicExtrapolate(gammaPrev,JPrev,slopePrev,gammaCur,JCur,slopeCur) ;
    gammaNew=max(gammaNew,xtrapf*gammaCur) ;
    gammaNew=min(gammaNew,gammaMax) ;

    gammaPrev=gammaCur ; JPrev=JCur ; slopePrev=slopeCur ;
    gammaCur=gammaNew ;
    [JCur,slopeCur]=Evaluate(gammaCur) ;

end

%%

info.Armijo    = Jmin <= J0+c1*gammaMin*slope0 ;
info.Curvature = abs(info.slopeMin) <= -c2*slope0 ;

if info.InfoLevel>=1
    fprintf(' LineSearchWolfe: gamma=%-12.6g J=%-12.6g slope=%-12.6g nFunc=%-3i nGrad=%-3i  (%s) \n',...
        gammaMin,Jmin,info.slopeMin,info.nFuncEvaluations,info.nGradEval,info.Message)
end

%% nested functions

    function [J,slope]=Evaluate(gamma)
        [J,slope]=Func(gamma) ;
        CheckScalar(J,'the first output of Func') ;
        CheckScalar(slope,'the second output of Func') ;
        info.nFuncEvaluations=info.nFuncEvaluations+1 ;
        info.nGradEval=info.nGradEval+1 ;
        if info.InfoLevel>=10
            fprintf('   eval %-3i : gamma=%-12.6g J=%-12.6g slope=%-12.6g \n',...
                info.nFuncEvaluations,gamma,J,slope)
        end
    end

    function UpdateBest(gamma,J,slope)
        if ~isempty(J) && isfinite(J) && J<JBest
            gammaBest=gamma ; JBest=J ;
            if isempty(slope) ; slopeBest=nan ; else ; slopeBest=slope ; end
        end
    end

    function [gamma,J]=GiveUp
        gamma=gammaBest ; J=JBest ; info.slopeMin=slopeBest ;
        if gammaBest==0
            info.iExit=-2 ;
            info.Message='no point better than phi(0) found.' ;
        else
            info.iExit=2 ;
            info.Message='evaluation limit reached, returning best point found.' ;
        end
    end

    function [gamma,J]=Zoom(gLo,JLo,slopeLo,gHi,JHi)

        % [gLo,gHi] brackets a minimiser. gLo has the lower J of the two and a
        % slope pointing into the bracket; gHi need not be the larger of the two.
        % slopeLo is known, the slope at gHi is not tracked.

        slopeHi=nan ;

        while true

            if info.nFuncEvaluations >= info.MaxFuncEvaluations
                [gamma,J]=GiveUp ; return
            end

            gLow=min(gLo,gHi) ; gHigh=max(gLo,gHi) ;
            Width=gHigh-gLow ;

            if Width <= info.xtol*max(1,gHigh)
                gamma=gammaBest ; J=JBest ; info.slopeMin=slopeBest ;
                info.iExit=3 ;
                info.Message='bracket collapsed before the curvature condition was met.' ;
                return
            end

            % safeguarded interpolation, kept away from the bracket ends
            if isfinite(slopeHi)
                gamma=CubicInterpolate(gLo,JLo,slopeLo,gHi,JHi,slopeHi) ;
            else
                gamma=QuadraticInterpolate(gLo,JLo,slopeLo,gHi,JHi) ;
            end
            if ~isfinite(gamma) || gamma<=gLow || gamma>=gHigh
                gamma=gLow+0.5*Width ;
            end
            gamma=min(max(gamma,gLow+0.1*Width),gHigh-0.1*Width) ;

            [J,slope]=Evaluate(gamma) ;
            UpdateBest(gamma,J,slope) ;

            if J > J0+c1*gamma*slope0 || J >= JLo
                gHi=gamma ; JHi=J ; slopeHi=slope ;
            else
                if abs(slope) <= -c2*slope0
                    info.slopeMin=slope ; info.iExit=1 ;
                    info.Message='strong Wolfe conditions satisfied.' ;
                    return
                end
                if slope*(gHi-gLo) >= 0
                    gHi=gLo ; JHi=JLo ; slopeHi=slopeLo ;
                end
                gLo=gamma ; JLo=J ; slopeLo=slope ;
            end

        end

    end

end

%% subfunctions

function CheckScalar(x,Name)

% The slope wanted here is the directional derivative dJ/dgamma, a scalar, and
% not the parameter gradient dJ/dp. If a gradient vector is passed in by mistake
% every subsequent if statement silently becomes an all() over its components,
% so this is caught up front.

if ~isnumeric(x) || ~isscalar(x)
    error('LineSearchWolfe:NotAScalar',...
        ['%s must be a numeric scalar, but has size [%s].\n',...
        'If this is a gradient with respect to the parameters, contract it with the\n',...
        'search direction first: phi''=g''*d for the l2 gradient, or phi''=dJdp''*(G*d)\n',...
        'for a Riesz-mapped gradient. See the header of LineSearchWolfe.m .'],...
        Name,strtrim(sprintf('%i ',size(x))))
end

if ~isfinite(x)
    error('LineSearchWolfe:NotFinite','%s must be finite, but is %g.',Name,x)
end

end

%%

function info=SetDefault(info,Field,Value)

if ~isfield(info,Field) || isempty(info.(Field))
    info.(Field)=Value ;
end

end

%%

function gamma=CubicInterpolate(a,fa,ga,b,fb,gb)

% minimiser of the cubic through (a,fa,ga) and (b,fb,gb), Nocedal & Wright 3.59

d1=ga+gb-3*(fa-fb)/(a-b) ;
Rad=d1*d1-ga*gb ;
if Rad<0 ; gamma=nan ; return ; end
d2=sign(b-a)*sqrt(Rad) ;
Den=gb-ga+2*d2 ;
if Den==0 ; gamma=nan ; return ; end
gamma=b-(b-a)*(gb+d2-d1)/Den ;

end

%%

function gamma=QuadraticInterpolate(a,fa,ga,b,fb)

% minimiser of the quadratic through (a,fa,ga) and (b,fb)

Den=fb-fa-ga*(b-a) ;
if Den<=0 ; gamma=nan ; return ; end
gamma=a-ga*(b-a)^2/(2*Den) ;

end

%%

function gamma=CubicExtrapolate(a,fa,ga,b,fb,gb)

% cubic extrapolation beyond b, used in the bracketing phase

gamma=CubicInterpolate(a,fa,ga,b,fb,gb) ;
if ~isfinite(gamma) || gamma<=b
    gamma=2*b ;
end

end
