function [gamma,f,Out,Info]=NewtonLineSearch(Func,f0,Out0,slope0,CtrlVar)

%% Line search along the Newton direction: backtracking with the Armijo condition
%
%   [gamma,f,Out,Info]=NewtonLineSearch(Func,f0,Out0,slope0,CtrlVar)
%
% (15 Oct 2026) Replaces rLineminUa.m (uv solver) and rLineminUa2.m (uvh solver).
%
% Func    : function handle, [f,Out]=Func(gamma), returning the merit function f at the step length gamma, and a structure Out with
%           any further quantities that are needed at the accepted step (eg rForce, rWork, D2, rBlocks). Func is evaluated at most
%           once for each step length, and the returned Out always corresponds to the returned gamma, so that no further evaluation is
%           required by the caller.
% f0,Out0 : the merit function and the further quantities at gamma=0. These are provided by the caller, as they can be calculated
%           from the residual already assembled for the Newton system, without a further evaluation of Func.
% slope0  : the slope of f at gamma=0. For the force residual (a weighted sum of squares of the residuals) and the Newton direction
%           this is slope0=-2*f0. If slope0 is not negative (or not finite), only actual reductions of f are accepted, and the step is
%           halved in each backtracking step.
%
% The full Newton step, gamma=1, is tried first. It is accepted if it fulfils the Armijo condition
%
%     f(gamma) <= f0 + c gamma slope0 ,      c=CtrlVar.NewtonLineSearch.ArmijoParameter (default 0.25).
%
% With slope0=-2*f0, c=0.25 means that the full step is accepted if it at least halves f, while for shorter steps the required
% reduction decreases in proportion to gamma. (The textbook value c=1e-4 accepts full steps that reduce f only marginally, although a
% shorter step often reduces f much more, and in tests resulted in considerably more Newton-Raphson iterations.) For c<0.5 full Newton
% steps are always accepted in the region of quadratic convergence.
%
% Otherwise the step length is reduced (backtracking). The new step length is the minimum of the quadratic through f0, slope0 and
% f(gamma) in the first backtracking step, and of the cubic through f0, slope0 and the last two function values in later steps (eg
% Nocedal and Wright, Numerical Optimization, 2006, Section 3.5). It is safeguarded to the interval [0.1,0.5]*gamma. Backtracking ends
% when the Armijo condition is fulfilled, when the next step length would be smaller than CtrlVar.NewtonLineSearch.GammaMin (default
% 1e-3), or after CtrlVar.NewtonLineSearch.MaxBacktrackingSteps (default 15) steps.
%
% Optionally (CtrlVar.NewtonLineSearch.Extrapolation=true, default false), step lengths larger than one are tried if the full step is
% accepted and the quadratic model through f0, slope0 and f(1) has its minimum beyond gamma=1. The step length is then increased by
% CtrlVar.NewtonLineSearch.ExtrapolationFactor (default 2) for as long as f decreases, at most MaxExtrapolationSteps (default 3) times.
%
% Refinement (only after backtracking): once the Armijo condition is fulfilled after backtracking, the minimum of the parabola through
% (0,f0) and the last two evaluated points is tried, and the better step is kept. This is repeated as long as each refinement reduces f
% by at least a factor of CtrlVar.NewtonLineSearch.RefinementRatio (default 0.5), at most CtrlVar.NewtonLineSearch.RefinementSteps
% (default 2) times. This brings the step closer to the minimum along the Newton direction where the full step was poor, and reduces the
% number of Newton-Raphson iterations. If the full step is accepted directly, no refinement is done.
%
% The best step found (smallest f) is returned. If no step reduces f, gamma=0, f=f0 and Out=Out0 are returned.
%
% Info.Converged : true if the Armijo condition is fulfilled, or if at least a clear reduction of f (by more than a fraction of 1e-5)
%                  was found with a step length of at least GammaMin. Otherwise false.
% Info.iarm      : number of backtracking steps
% Info.nFuncEval : number of evaluations of Func
% Info.nExtrapolationSteps : number of extrapolation steps
% Info.nRefinementSteps : number of refinement steps
% Info.Infovector: [gamma f] of all evaluated step lengths, including gamma=0, sorted by gamma
% Info.Direction : "N " (Newton direction), for compatibility with the output of the Newton-Raphson solvers
%
%%

narginchk(5,5)

%% options
c=0.25 ; gammaMin=1e-3 ; nBacktrackMax=15 ; Extrapolation=false ; nExtrapolationMax=3 ; ExtrapolationFactor=2 ; InfoLevel=0 ;
nRefineMax=2 ; RefineRatio=0.5 ;
if isfield(CtrlVar,"NewtonLineSearch")
    Opt=CtrlVar.NewtonLineSearch ;
    if isfield(Opt,"ArmijoParameter")       ; c=Opt.ArmijoParameter ; end
    if isfield(Opt,"GammaMin")              ; gammaMin=Opt.GammaMin ; end
    if isfield(Opt,"MaxBacktrackingSteps")  ; nBacktrackMax=Opt.MaxBacktrackingSteps ; end
    if isfield(Opt,"Extrapolation")         ; Extrapolation=Opt.Extrapolation ; end
    if isfield(Opt,"MaxExtrapolationSteps") ; nExtrapolationMax=Opt.MaxExtrapolationSteps ; end
    if isfield(Opt,"ExtrapolationFactor")   ; ExtrapolationFactor=Opt.ExtrapolationFactor ; end
    if isfield(Opt,"RefinementSteps")       ; nRefineMax=Opt.RefinementSteps ; end
    if isfield(Opt,"RefinementRatio")       ; RefineRatio=Opt.RefinementRatio ; end
end
if isfield(CtrlVar,"InfoLevelBackTrack") ; InfoLevel=CtrlVar.InfoLevelBackTrack ; end

isSlope = isfinite(slope0) && slope0<0 ;

%% evaluated points (gamma=0 is included, with the values provided by the caller)
gList=0 ; fList=f0 ; OutList={Out0} ;
nEval=0 ; iarm=0 ; nExtrapolation=0 ; nRefine=0 ; ArmijoFulfilled=false ;

Info.Converged=false ; Info.iarm=0 ; Info.nFuncEval=0 ; Info.nExtrapolationSteps=0 ; Info.nRefinementSteps=0 ; Info.Infovector=[0 f0] ; Info.Direction="N " ;

if ~(isfinite(f0) && f0>0)   % nothing to be reduced (f0=0), or f0 is not usable
    gamma=0 ; f=f0 ; Out=Out0 ; Info.Converged=(f0==0) ;
    return
end

%% full Newton step
g=1 ;
[fg,Og]=EvalFunc(g) ;

%% backtracking
gPrev=NaN ; fPrev=NaN ;
while true

    if Armijo(g,fg)
        ArmijoFulfilled=true ;
        break
    end

    if iarm>=nBacktrackMax
        break
    end

    gNew=NewStepLength(g,fg,gPrev,fPrev) ;

    if gNew<gammaMin
        break
    end

    gPrev=g ; fPrev=fg ;
    g=gNew ; iarm=iarm+1 ;
    [fg,Og]=EvalFunc(g) ;

end

%% refinement, only after backtracking
if ArmijoFulfilled && iarm>0 && nRefineMax>0
    gA=g ; fA=fg ; gB=gPrev ; fB=fPrev ;     % the accepted step, and the previous (larger) step
    while nRefine<nRefineMax && isfinite(fA) && isfinite(fB)
        gR=ParabolaMin(gA,fA,gB,fB) ;      % minimum of the parabola through (0,f0), (gA,fA) and (gB,fB)
        if ~(isfinite(gR) && gR>0) || abs(gR-gA)<0.1*gA
            break
        end
        [fR,~]=EvalFunc(gR) ;
        nRefine=nRefine+1 ;
        if fR<fA
            ImprovedEnough = fR < RefineRatio*fA ;
            if gR<gA ; gB=gA ; fB=fA ; end % keep a bracket around the minimum
            gA=gR ; fA=fR ;
            if ~ImprovedEnough ; break ; end
        else
            break
        end
    end
end

%% optional extrapolation, only if the full Newton step was accepted
if Extrapolation && ArmijoFulfilled && g==1 && isSlope
    gq=-slope0/(2*(fg-f0-slope0)) ;     % minimum of the quadratic through f0, slope0 and f(1)
    if isfinite(gq) && gq>1 && fg-f0-slope0>0
        fBest=fg ;
        while nExtrapolation<nExtrapolationMax
            gE=g*ExtrapolationFactor ;
            [fE,~]=EvalFunc(gE) ;
            nExtrapolation=nExtrapolation+1 ;
            if fE<fBest
                g=gE ; fBest=fE ;
            else
                break
            end
        end
    end
end

%% the best step found
[fBest,I]=min(fList) ;
gamma=gList(I) ; f=fBest ; Out=OutList{I} ;

if gamma==0
    Converged=false ;
else
    Converged = Armijo(gamma,f) || ( (f0-f)/f0 > 1e-5 && gamma>=gammaMin ) ;
end

[~,J]=sort(gList) ;
Info.Converged=Converged ;
Info.iarm=iarm ;
Info.nFuncEval=nEval ;
Info.nExtrapolationSteps=nExtrapolation ;
Info.nRefinementSteps=nRefine ;
Info.Infovector=[gList(J)' fList(J)'] ;

if InfoLevel>=10
    fprintf(" NewtonLineSearch: gamma=%-10.5g f/f0=%-10.5g | Armijo fulfilled=%i | backtracking steps=%i | refinement steps=%i | extrapolation steps=%i | evaluations=%i \n",...
        gamma,f/f0,ArmijoFulfilled,iarm,nRefine,nExtrapolation,nEval) ;
end


%% nested functions

    function [fv,Ov]=EvalFunc(gv)
        % evaluates Func at gv, and keeps the result. Non-finite values are treated as a failure (f=inf).
        [fv,Ov]=Func(gv) ;
        fv=full(fv) ;
        if ~isfinite(fv) || ~isreal(fv)
            fv=inf ;
        end
        nEval=nEval+1 ;
        gList(end+1)=gv ; fList(end+1)=fv ; OutList{end+1}=Ov ;
    end

    function isOK=Armijo(gv,fv)
        if isSlope
            isOK = fv <= f0+c*gv*slope0 ;
        else
            isOK = fv < f0 ;
        end
    end

    function gm=ParabolaMin(g1,f1,g2,f2)
        % minimum of the parabola through (0,f0), (g1,f1) and (g2,f2). NaN if the parabola has no minimum.
        gm=NaN ;
        den=g1*g2*(g1-g2) ;
        if den==0 ; return ; end
        aa=( g2*(f1-f0)-g1*(f2-f0))/den ;          % coefficient of gamma^2
        bb=(-g2^2*(f1-f0)+g1^2*(f2-f0))/den ;       % coefficient of gamma
        if aa>0
            gm=-bb/(2*aa) ;
        end
    end

    function gn=NewStepLength(gv,fv,gp,fp)
        % new step length from quadratic (first step) or cubic (later steps) interpolation, safeguarded to [0.1,0.5]*gv
        gn=NaN ;
        if isSlope && isfinite(fv)
            if isnan(gp) || ~isfinite(fp)
                % quadratic through f0, slope0 and f(gv)
                den=2*(fv-f0-slope0*gv) ;
                if den>0
                    gn=-slope0*gv^2/den ;
                end
            else
                % cubic through f0, slope0, f(gv) and f(gp)
                d1=fv-f0-slope0*gv ;
                d2=fp-f0-slope0*gp ;
                den=gv^2*gp^2*(gv-gp) ;
                a=( gp^2*d1-gv^2*d2)/den ;
                b=(-gp^3*d1+gv^3*d2)/den ;
                if a==0
                    if b~=0 ; gn=-slope0/(2*b) ; end
                else
                    disc=b^2-3*a*slope0 ;
                    if disc>=0
                        gn=(-b+sqrt(disc))/(3*a) ;
                    end
                end
            end
        end
        if ~isfinite(gn) || ~isreal(gn)
            gn=0.5*gv ;
        end
        gn=min(max(gn,0.1*gv),0.5*gv) ;
    end

end