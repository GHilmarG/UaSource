function [gammamin,rmin,du,dv,dh,dl,BackTrackInfo,rForce,rWork,D2,rBlocks] = rLineminUa2(CtrlVar,UserVar,func,r0,r1,K,L,du0,dv0,dh0,dl0,dJdu,dJdv,dJdh,dJdl,Normalisation,M)

%%
%
% Line search for the non-linear (Newton) iteration of the uv and uvh problems.
%
%   r(gamma) = func(gamma,Du,Dv,Dh,Dl)
%
% is the normalised squared norm of the KKT residual at the point displaced by gamma*[Du;Dv;Dh;Dl].
%
% The Newton direction is the solution of
%
%   [ K  L' ]  [dx0]  = [ dJdx ]   =: R
%   [ L  0  ]  [dl0]    [ dJdl ]
%
% so that, with H=[K L'; L 0], the Gauss-Newton model of r along a direction s is  Q(gamma)=sum(w.*(R-gamma*H*s).^2),
% where w holds the weight of each row of the residual in the merit function r. The weights are 1/Normalisation for all rows if the
% residuals are normalised together, and block-dependent if CtrlVar.uvhResidualNormalisation="blockwise". They are obtained from the block
% residuals returned by func (see ResidualWeights below), so that the model is consistent with r. (Only pooled weights for the uv problem.)
%
% The search proceeds in three stages. Each stage is only entered if the previous ones have not given a satisfactory reduction.
%
%  1) Newton:   the full Newton step with backtracking (min step size CtrlVar.BacktrackingGammaMin=0.001).
%               The result is accepted unless the Newton step failed (no reduction in r), or the reduction is weak, ie
%               r/r0 > rRatioWeak (default 0.9).
%
%  2) Cauchy:   (cheap, no linear solve with K) steepest-descent-type step  s = Dm^(-1) R, projected on the constraints, where
%               Dm is a metric matrix. Two metrics are tried: Dm=diag(|K|), ie the diagonal of the Jacobian, which has the correct units for
%               u, v and h ("D"), and the block mass matrix ("M"). CtrlVar.rLineMinUaCauchy = "DM" (default, both are tried, the better is used), "D", "M", or "none".
%               The step length is the minimiser of the Gauss-Newton model along s, followed by backtracking on the true r.
%               Accepted if r/r0 < rRatioGood (default 0.7).
%
%  3) LM:       (OFF by default) Levenberg-Marquardt type step   (K + lambda*diag|K|) s = R   (same KKT structure as the Newton step).
%               lambda=0 is the Newton step and large lambda tends to a diagonally scaled steepest descent step. Starting with the
%               value of lambda that worked last time (persistent), lambda is increased by a factor of 10 (at most 3 trials, lambda<=100) until a step with
%               r/r0 < rRatioGood is found, stopping if a trial is worse than the previous one. Each trial requires one linear solve. lambda is then reduced by a factor 0.3
%               for the next call if the full LM step was accepted. Switch on with CtrlVar.rLineMinUaLM=true. In tests on a Greenland case it gave no improvement over
%               the Cauchy steps, at the cost of extra solves. Optional: CtrlVar.rLineMinUaLMTrigger (try LM if best r/r0 is above this value), CtrlVar.rLineMinUaLMScaling ("diag" or "rowsum").
%
% The best step found is returned. If no step reduces r then gammamin=0 and the returned step is zero.
%
% Returned BackTrackInfo.Direction is "N " (Newton), "CD" or "MD" (Cauchy with D or M metric) or "LM". Options: CtrlVar.rLineMinUaWeakRatio, rLineMinUaGoodRatio.
%
% The meaning of the returned gammamin: Newton: the step length. Cauchy: step length divided by the model minimiser.
%                                       LM: step length with respect to the LM step.
%%

persistent lambdaMemory
if isempty(lambdaMemory) ; lambdaMemory=1 ; end

%% options
rRatioWeak=0.9 ; rRatioGood=0.7 ; nLMmax=3 ; lambdaUp=10 ; lambdaDown=0.3 ; lambdaMin=1e-4 ; lambdaMax=1e8 ;
CauchyMetric="DM" ; LMOn=false ; LMScaling="diag" ; lambdaCap=100 ;
if isfield(CtrlVar,"rLineMinUaWeakRatio")  ; rRatioWeak=CtrlVar.rLineMinUaWeakRatio ; end
if isfield(CtrlVar,"rLineMinUaGoodRatio")  ; rRatioGood=CtrlVar.rLineMinUaGoodRatio ; end
if isfield(CtrlVar,"rLineMinUaCauchy")     ; CauchyMetric=string(CtrlVar.rLineMinUaCauchy) ; end
if isfield(CtrlVar,"rLineMinUaLM")         ; LMOn=logical(CtrlVar.rLineMinUaLM) ; end
if isfield(CtrlVar,"rLineMinUaLMmaxTrials"); nLMmax=CtrlVar.rLineMinUaLMmaxTrials ; end
if isfield(CtrlVar,"rLineMinUaLMScaling")   ; LMScaling=string(CtrlVar.rLineMinUaLMScaling) ; end
rRatioLMTrigger=rRatioGood ; if isfield(CtrlVar,"rLineMinUaLMTrigger") ; rRatioLMTrigger=CtrlVar.rLineMinUaLMTrigger ; end
if isfield(CtrlVar,"rLineMinUaLMmaxLambda") ; lambdaCap=CtrlVar.rLineMinUaLMmaxLambda ; end

%% these will be the returns if no reduction is found
gammamin=0 ; rmin=r0 ;  du=du0*0 ; dv=dv0*0 ; dh=dh0*0 ; dl=dl0*0 ;
rForce=r0 ; rWork=nan ; D2=nan ; rBlocks=[nan nan nan nan] ;

%% variables
if isempty(dh0)
    nBlk=2 ; Rx=[dJdu;dJdv] ; sol0=[du0;dv0] ;
else
    nBlk=3 ; Rx=[dJdu;dJdv;dJdh] ; sol0=[du0;dv0;dh0] ;
end
nM=size(M,1) ; nx=numel(Rx) ; nL=size(L,1) ;
dl0=dl0(:) ; dJdl=dJdl(:) ;
Rfull=[Rx;dJdl] ;
NewtonDir=[sol0;dl0] ;

dbg=isfield(CtrlVar,"rLineMinUaDebug") && CtrlVar.rLineMinUaDebug ;
EvalCounter("reset") ; nSolve=0 ;
funcC=@(varargin) CountedCall(func,varargin{:}) ;
Rfun=@(gamma,s) Eval(funcC,gamma,s,nM,nBlk) ;      % scalar r at displacement gamma*s

%% 1) Newton step with backtracking
rNewtonFunc=@(gamma) Rfun(gamma,NewtonDir) ;
if isnan(r0) || isempty(r0) ; r0=rNewtonFunc(0) ; rmin=r0 ; end
if isnan(r1) || isempty(r1) ; r1=rNewtonFunc(1) ; end
CtrlVarN=CtrlVar ; CtrlVarN.BacktrackingGammaMin=0.001 ; CtrlVarN.BacktracFigName="Line Search in Newton Direction" ;
[gN,rN,BTN]=BackTracking(-2*r0,1,r0,r1,rNewtonFunc,CtrlVarN) ;

best.r=r0 ; best.s=NewtonDir*0 ; best.gamma=0 ; best.label="  " ; best.info=BTN ; best.lambda=nan ;
if rN < best.r
    best.r=rN ; best.s=gN*NewtonDir ; best.gamma=gN ; best.label="N " ; best.info=BTN ;
end

NewtonFailed = isnan(rN) || ~(rN<r0) || ~BTN.Converged ;
TryFallback  = NewtonFailed || rN/r0 > rRatioWeak || (LMOn && rN/r0 > rRatioLMTrigger) ;
rCauchy=nan ; rLM=nan ; lambdaUsed=nan ;

if TryFallback

    if nL>0 ; H=[K L.' ; L sparse(nL,nL)] ; else ; H=K ; end
    dK=abs(full(diag(K))) ; dK=max(dK,1e-12*max(dK)) ; if ~all(isfinite(dK)) || max(dK)==0 ; dK=ones(nx,1) ; end

    %% Weights of the rows of the residual in the merit function (consistent with the way r is calculated in func)
    [wRow,wMode]=ResidualWeights(funcC,NewtonDir,nM,nBlk,nx,nL,Rfull,r0,Normalisation) ;
    if dbg ; fprintf("   [weights] %s, w_u=%.3g w_v=%.3g w_h=%.3g (1/Normalisation=%.3g)\n",wMode,wRow(1),wRow(nM+1),wRow(min(2*nM+1,numel(wRow))),1/Normalisation) ; end

    %% 2) Cauchy step (cheap)
    if CauchyMetric=="DM" ; metricList=["D" "M"] ; elseif CauchyMetric=="none" ; metricList=strings(1,0) ; else ; metricList=CauchyMetric ; end
    for metric=metricList
        if best.r/r0 < rRatioGood ; break ; end

        if metric=="D"
            if nL==0
                sx=Rx./dK ; dlC=zeros(0,1) ;
            else
                [sx,dlC]=solveKApeSymmetric(spdiags(dK,0,nx,nx),L,Rx,dJdl,sol0,dl0,CtrlVar) ;
            end
            labelC="CD" ;
        else
            Mb=kron(speye(nBlk),M) ;
            [sx,dlC]=solveKApeSymmetric(Mb,L,Rx,dJdl,sol0,dl0,CtrlVar) ;
            labelC="MD" ;
        end
        sC=[sx(:);dlC(:)] ;
        HsC=H*sC ; num=Rfull'*(wRow.*HsC) ; den=HsC'*(wRow.*HsC) ;

        if num>0 && den>0
            gC=num/den ; slopeC=-2*num ;      % r(gamma) already contains the weights
            rC1=Rfun(gC,sC) ;
            if dbg   % compare the model slope with a finite-difference slope of the true r, and with the slope based on pooled weights
                gFD=0.01*gC ; slopeFD=(Rfun(gFD,sC)-r0)/gFD ; slopePooled=-2*(Rfull'*HsC)/Normalisation ;
                fprintf("   [Cauchy %s slope check] model slope=%.4g | finite-difference slope=%.4g | slope with pooled Normalisation=%.4g | model/FD=%.4g pooled/FD=%.4g\n",labelC,slopeC,slopeFD,slopePooled,slopeC/slopeFD,slopePooled/slopeFD) ;
            end
            CtrlVarC=CtrlVar ; CtrlVarC.NewtonAcceptRatio=0.9 ; CtrlVarC.BacktrackingGammaMin=1e-10 ; CtrlVarC.LineSearchAllowedToUseExtrapolation=false ;
            CtrlVarC.BacktracFigName="Line Search in Cauchy direction" ;
            [gCm,rCtrial,BTC]=BackTracking(slopeC,gC,r0,rC1,@(g) Rfun(g,sC),CtrlVarC) ;
            if dbg ; fprintf("   [Cauchy %s] gC=%.3g  r(gC)/r0=%.4g -> backtracked g=%.3g r/r0=%.4g  (|sC|/|N|=%.3g)\n",labelC,gC,rC1/r0,gCm,rCtrial/r0,norm(sC)/norm(NewtonDir)) ; end
            if isnan(rCauchy) || rCtrial<rCauchy ; rCauchy=rCtrial ; end
            if rCtrial < best.r
                best.r=rCtrial ; best.s=gCm*sC ; best.gamma=gCm/gC ; best.label=labelC ; best.info=BTC ;
            end
        end
    end

    %% 3) Levenberg-Marquardt type step
    if LMOn && best.r/r0 >= rRatioLMTrigger

        if LMScaling=="rowsum"
            dLM=full(sum(abs(K),2)) ; dLM=max(dLM,1e-12*max(dLM)) ;     % abs row sums: K+lambda*D is strictly diagonally dominant for lambda>1
        else
            dLM=dK ;
        end
        DD=spdiags(dLM,0,nx,nx) ;
        lambda=min(max(lambdaMemory,lambdaMin),lambdaCap) ; success=false ; fullStep=false ; rPrev=inf ; rBestLM=inf ; lambdaBestLM=lambda ;
        CtrlVarL=CtrlVar ; CtrlVarL.NewtonAcceptRatio=min(rRatioGood,rRatioLMTrigger) ; CtrlVarL.BacktrackingGammaMin=1e-3 ; CtrlVarL.LineSearchAllowedToUseExtrapolation=false ;
        CtrlVarL.BacktracFigName="Line Search in LM direction" ;

        for iLM=1:nLMmax

            [sx,dlL]=solveKApe(K+lambda*DD,L,Rx,dJdl,sol0,dl0,CtrlVar) ; nSolve=nSolve+1 ;
            sL=[sx(:);dlL(:)] ;
            HsL=Rfull ; HsL(1:nx)=HsL(1:nx)-lambda*(dLM.*sx(:)) ;     % H*sL, using (K+lambda D) sx + L'dl = Rx
            slopeL=-2*(Rfull'*(wRow.*HsL)) ;
            rLtrial=inf ;

            if slopeL<0 && isfinite(slopeL)
                rL1=Rfun(1,sL) ;
                [gL,rLtrial,BTL]=BackTracking(slopeL,1,r0,rL1,@(g) Rfun(g,sL),CtrlVarL) ;
                if dbg ; fprintf("   [LM trial %d] lambda=%.3g  slope0=%.3g  |sL|/|N|=%.3g  r(1)/r0=%.4g -> backtracked g=%.3g r/r0=%.4g\n",iLM,lambda,slopeL,norm(sL)/norm(NewtonDir),rL1/r0,gL,rLtrial/r0) ; end
                if rLtrial < best.r
                    best.r=rLtrial ; best.s=gL*sL ; best.gamma=gL ; best.label="LM" ; best.info=BTL ; best.lambda=lambda ;
                end
                if rLtrial < rBestLM ; rBestLM=rLtrial ; lambdaBestLM=lambda ; rLM=rLtrial ; lambdaUsed=lambda ; end
                if rLtrial/r0 < min(rRatioGood,rRatioLMTrigger)
                    success=true ; fullStep=(gL>0.99) ; break
                end
            elseif dbg
                fprintf("   [LM trial %d] lambda=%.3g  not a descent direction (slope0=%.3g)\n",iLM,lambda,slopeL) ;
            end

            % Stop going up in lambda if this trial is not better than the previous one. Larger lambda gives smaller and smaller steps.
            if iLM>1 && rLtrial>=rPrev ; break ; end
            rPrev=rLtrial ;
            lambda=lambda*lambdaUp ;
            if lambda>lambdaCap ; break ; end
        end

        if success
            if fullStep ; lambdaMemory=max(lambda*lambdaDown,lambdaMin) ; else ; lambdaMemory=lambda ; end
        elseif isfinite(rBestLM) && rBestLM<r0
            lambdaMemory=lambdaBestLM ;
        else
            lambdaMemory=3 ;
        end
    end
end

%% return values
NoReduction=~(best.r<r0) ;
if ~NoReduction
    s=best.s ;
    du=s(1:nM) ; dv=s(nM+1:2*nM) ;
    if nBlk==3 ; dh=s(2*nM+1:3*nM) ; end
    dl=s(nBlk*nM+1:end) ;
    gammamin=best.gamma ; rmin=best.r ;
    try
        [~,~,~,rForce,rWork,D2,rBlocks]=Eval(funcC,1,s,nM,nBlk) ;
    catch
        [~,~,~,rForce,rWork,D2]=Eval(funcC,1,s,nM,nBlk) ;
    end
    if best.label=="CD" || best.label=="MD" ; rWork=nan ; D2=nan ; end    % no meaning for these steps
    BackTrackInfo=best.info ; BackTrackInfo.Direction=best.label ;
else
    BackTrackInfo=best.info ;
    BackTrackInfo.Direction="  " ;
end

rRatioMin=0.99999 ;
if NoReduction || rmin/r0 > rRatioMin
    BackTrackInfo.Converged=false ;
else
    BackTrackInfo.Converged=true ;
end

if CtrlVar.InfoLevelNonLinIt >= 5
    fprintf(" rLineminUa: dir=%s  r/r0=%-10.4g gamma=%-10.4g | Newton r/r0=%-10.4g | Cauchy r/r0=%-10.4g | LM r/r0=%-10.4g lambda=%-8.3g | #eval=%i #solve=%i \n",...
        best.label,rmin/r0,gammamin,rN/r0,rCauchy/r0,rLM/r0,lambdaUsed,EvalCounter("get"),nSolve) ;
end

end


function varargout=Eval(func,gamma,s,nM,nBlk)
% evaluates func at the displacement gamma*s, where s=[du;dv;(dh);dl]
if nBlk==3
    [varargout{1:max(nargout,1)}]=func(gamma,s(1:nM),s(nM+1:2*nM),s(2*nM+1:3*nM),s(3*nM+1:end)) ;
else
    [varargout{1:max(nargout,1)}]=func(gamma,s(1:nM),s(nM+1:2*nM),s(2*nM+1:end)) ;
end
end

function varargout=CountedCall(f,varargin)
EvalCounter("inc") ;
[varargout{1:max(nargout,1)}]=f(varargin{:}) ;
end

function n=EvalCounter(cmd)
persistent count
if isempty(count) ; count=0 ; end
switch cmd
    case "reset" ; count=0 ;
    case "inc"   ; count=count+1 ;
end
n=count ;
end


function [w,wMode]=ResidualWeights(func,NewtonDir,nM,nBlk,nx,nL,Rfull,r0,Normalisation)

%%
% Weight of each row of the residual vector in the merit function r. 
%
%   r = sum_b  w_b |R_b|^2,    b=u,v,h,l
%
% func returns, as its seventh output, the block residuals rBlocks_b = w_b |R_b|^2 (for the uvh problem), so w_b=rBlocks_b/|R_b|^2.
% If this can not be done, or if the weighted sum of squares does not reproduce r0 (for example if r is the work residual), then
% the weights revert to the pooled value 1/Normalisation.
%%

w=zeros(nx+nL,1)+1/Normalisation ;
wMode="pooled" ;

if nBlk~=3 ; return ; end       % the uv problem only has pooled residuals (and its 7th output is not rBlocks)

try
    [~,~,~,~,~,~,rB]=Eval(func,0,NewtonDir,nM,nBlk) ;
catch
    return
end

if ~isnumeric(rB) || numel(rB)~=4 ; return ; end
rB=full(rB(:)) ;

rows={1:nM , nM+1:2*nM , 2*nM+1:3*nM , 3*nM+1:nx+nL} ;
wb=nan(4,1) ;
for b=1:4
    R2=sum(Rfull(rows{b}).^2) ;
    if R2>0 && rB(b)>0 ; wb(b)=rB(b)/R2 ; end
end
% u and v share their weights, and a block with zero residual does not contribute (its weight is irrelevant if the step satisfies the constraints)
if isnan(wb(1)) ; wb(1)=wb(2) ; end
if isnan(wb(2)) ; wb(2)=wb(1) ; end
wb(isnan(wb))=1/Normalisation ;

wTest=zeros(nx+nL,1) ;
for b=1:4 ; wTest(rows{b})=wb(b) ; end

rModel=sum(wTest.*Rfull.^2) ;
if abs(rModel-r0) <= 1e-6*abs(r0)
    w=wTest ;
    if max(wb(1:3))>1.0001*min(wb(1:3)) ; wMode="blockwise" ; else ; wMode="pooled" ; end
else
    wMode="pooled (fallback: weighted residuals do not reproduce r0)" ;
end

end
