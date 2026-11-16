function [p,UserVar,RunInfo]=UaOptimisationHessianVectorBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub,HessianSetup)

%%
%
%   [p,UserVar,RunInfo]=UaOptimisationHessianVectorBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub,HessianSetup)
%
% Inexact Gauss-Newton (default) or Newton trust-region method using Hessian-vector products, CtrlVar.Inverse.MinimisationMethod=
% "-UaOptimisation-HessianVectorBased-" (introduced 15 November 2026 as "-UaOptimisation-NewtonCG-", which is still accepted).
%
% The cost function (Gaussian errors and Gaussian priors) is a weighted nonlinear least-squares problem. With
% CtrlVar.Inverse.HVP.Type="-GN-" (default) the Gauss-Newton Hessian is used, ie the terms proportional to the residuals (the
% second derivatives of the forward model weighted by the adjoint variables) are omitted; the Gauss-Newton Hessian is positive
% semi-definite. With "-full-" the full (Newton) Hessian is used. The Hessian is never formed, only Hessian-vector products.
%
% At each (outer) iteration the Newton system H s = -g is solved approximately with the conjugate-gradient method, preconditioned
% with the metric G (the Hessian of the regularisation term, MUA.G), subject to the trust-region constraint ||s||_G <= Delta
% (Steihaug-Toint). The CG iteration stops when
%
%   (1) the preconditioned residual is small enough: sqrt(r' G^{-1} r) <= eta sqrt(g' G^{-1} g), eta from Eisenstat-Walker (choice 1)
%   (2) the step reaches the trust-region boundary,
%   (3) a direction of negative curvature is found (possible with the full Hessian), the step then goes to the boundary, or
%   (4) the maximum number of CG iterations is reached.
%
% The Hessian-vector products are provided by HessianSetup(x,F,JGHouts), see HessianVectorProductSetup.m. Each product costs a linearised
% forward and an adjoint solve. By default the Gauss-Newton Hessian is used (CtrlVar.Inverse.HVP.Type).
%
% The step is accepted if rho=(actual decrease)/(predicted decrease) > 1e-4, and the trust-region radius is reduced (rho<0.25) or
% increased (rho>0.75 and the step on the boundary) as usual.
%
% *Box constraints* (CtrlVar.Inverse.HessianVectorBased.Bounds):
%
%   "-projected-"    (default) : the bounds are imposed explicitly, in the physical parameters (no box transform). Each step consists
%                                of (1) a Cauchy step along the projected steepest-descent path P(x-t g), giving sufficient decrease of
%                                the quadratic model within the trust region, which identifies the variables at their bounds (the
%                                active set), (2) a Steihaug-Toint CG step in the remaining (free) variables, preconditioned with G
%                                restricted to the free variables, and (3) a projected search along that step, accepting it if the
%                                model decrease is at least that of the Cauchy step (cf. TRON, Lin and More, 1999). The quadratic model
%                                is then exactly the Gauss-Newton (or full) model in the physical parameters.
%   "-BoxTransform-"           : the bounds are imposed through the box transform (CtrlVar.Inverse.BoxTransform=true), so that the
%                                problem is unconstrained. The preconditioner G is used as it is for the transformed variables. Close to
%                                the bounds the transform is strongly nonlinear, which can limit the convergence.
%
% CtrlVarValidityCheck.m sets CtrlVar.Inverse.BoxTransform accordingly.
%
% *Nodal scaling of the trust region* (CtrlVar.Inverse.HessianVectorBased.NodalTrustRegionScaling=true, experimental): the trust region is
% ||s||_M <= Delta with M=G+diag(w). After a failed step, if the mismatch between the actual and the predicted change of the gradient
% is concentrated in a few variables, their weights are increased, so that the step is limited there but not elsewhere. Intended for
% problems where the quadratic model fails locally, for example close to flotation. M is then also the CG preconditioner.
%
% Settings: CtrlVar.Inverse.Iterations (maximum number of Newton iterations), CtrlVar.Inverse.HessianVectorBased.MaxInnerIterations,
% CtrlVar.Inverse.HessianVectorBased.InitialTrustRegionRadius (default sqrt(numel(p))), CtrlVar.Inverse.HessianVectorBased.MaxForcingTerm.
%
% see also: HessianVectorProductSetup.m, HessianVectorProduct.m, AssembleHessianBlocks.m, UaOptimisationGradientBased.m
%
%%

narginchk(9,9)

Settings=struct("MaxInnerIterations",200,"InitialTrustRegionRadius",[],"MaxForcingTerm",0.5,"Bounds","-projected-","NodalTrustRegionScaling",false,"LocalisationFraction",0.01,"LocalisationShare",0.5);
if isfield(CtrlVar.Inverse,"HessianVectorBased")
    fn=fieldnames(CtrlVar.Inverse.HessianVectorBased);
    for k=1:numel(fn) ; Settings.(fn{k})=CtrlVar.Inverse.HessianVectorBased.(fn{k}) ; end
end
UserGivenRadius=~isempty(Settings.InitialTrustRegionRadius) ;   % (16 Nov 2026) a given radius takes precedence over a saved one
if isempty(Settings.InitialTrustRegionRadius)
    Settings.InitialTrustRegionRadius=sqrt(numel(p));
end

G=MUA.G ;
if isfield(MUA,"dG") && ~isempty(MUA.dG)
    Prec=@(r) MUA.dG\r ;
else
    dG=decomposition(G) ;
    Prec=@(r) dG\r ;
end

isBinversion=contains(CtrlVar.Inverse.InvertFor,"-B-");
p=p(:);

% bounds
if isempty(plb) ; plb=-inf(size(p)) ; end
if isempty(pub) ; pub=+inf(size(p)) ; end
plb=plb(:) ; pub=pub(:) ;
Projected=any(isfinite(plb)) || any(isfinite(pub)) ;
if Projected && string(Settings.Bounds)~="-projected-"
    error("UaOptimisationHessianVectorBased:Bounds","Explicit bounds given, but CtrlVar.Inverse.HessianVectorBased.Bounds=%s. Use -projected-, or the box transform.",Settings.Bounds)
end
p=min(max(p,plb),pub) ;   % feasible starting point

%% initial point
[J,g,~,fOuts,Fx]=func(p);
g=g(:);
% (16 Nov 2026) Factorisations of G (or M) restricted to the free variables, kept while the set of free variables (and M) does not
% change. When all variables are free, the factorisation of G (MUA.dG) or of M is used directly. See RestrictedSolver below.
StepCache=struct("Free",[],"Key",[],"dF",[],"nFactorisations",0) ;   % for the CG step (G, or M with nodal weights)
NormCache=StepCache ;                                                       % for the (projected) gradient norm (G)
ModelNormCache=StepCache ;                                                  % for the norm of the gradient predicted by the model (G)

[gM,NormCache]=GradientNorm(p,g,G,Prec,plb,pub,Projected,NormCache);     % (projected) gradient norm in the dual (inverse G) norm
gM0=gM;
Delta=Settings.InitialTrustRegionRadius;
Delta0=Delta;

%% RunInfo
Fields=["Iterations","J","R","I","StepSize","GradNorm","Decrement","Juv","Jhdot","JBobs","InnerIterations","TrustRegionRadius"];
for ff=Fields
    if ~isfield(RunInfo.Inverse,ff) || isempty(RunInfo.Inverse.(ff))
        RunInfo.Inverse.(ff)=NaN;
    end
end
if isempty(RunInfo.Inverse.Iterations) || all(isnan(RunInfo.Inverse.Iterations)) || numel(RunInfo.Inverse.Iterations)<=1
    It0=0;
    RunInfo=StoreIteration(RunInfo,1,It0,J,fOuts,0,gM,NaN,NaN,NaN,Fields);
else
    It0=max(RunInfo.Inverse.Iterations);
    n=numel(RunInfo.Inverse.J);
    for ff=Fields % same length as J (fields introduced later, or not filled by other optimisers)
        if numel(RunInfo.Inverse.(ff))<n ; RunInfo.Inverse.(ff)=[RunInfo.Inverse.(ff)(:) ; nan(n-numel(RunInfo.Inverse.(ff)),1)] ; end
    end
end

fprintf("\n +++++++++++ At start of inversion (HessianVectorBased, %s Hessian):  J=%-g  I=%-g  R=%-g  |g|_G=%-g \n",HVPType(CtrlVar),J,fOuts.I,fOuts.R,gM)
fprintf("\n   It  inner  CG exit           J                I                R            ||g||_G       ||s||_G       Delta        rho")
if isBinversion ; fprintf("          JBobs") ; end
if Projected ; fprintf("   at bounds") ; end
if Settings.NodalTrustRegionScaling ; fprintf("  weighted") ; end
fprintf("      time\n")

nRejected=0;
tStart=tic;

% (15 November 2026) Nodal scaling of the trust region: ||s||_M <= Delta with M = G + diag(w). The weights w are zero to begin with,
% and are increased at the few variables where the quadratic model has failed (where the actual change in the gradient differs most
% from the change predicted by the model), see below. This keeps the steps small where the problem is strongly nonlinear (eg close
% to flotation), without limiting the step everywhere else.
n=numel(p) ; w=zeros(n,1) ;

% (15 November 2026) Forcing term (relative tolerance of the inner CG iteration), Eisenstat-Walker choice 1: eta measures how well
% the quadratic model predicted the gradient at the new point, eta = | |g_new| - |g_old + H s| | / |g_old|  (in the dual G norm),
% safeguarded and limited to MaxForcingTerm. This depends only on the last step, and not on the gradient at the start of the run,
% so that it also works for restart runs (previously eta=min(MaxForcingTerm,sqrt(|g|/|g0|)), with g0 the gradient at the start of
% the run, which after a restart close to the minimum gave far too loose inner tolerances).
eta=Settings.MaxForcingTerm ;

% (16 Nov 2026) Restart runs continue with the trust-region state at the end of the previous run (trust-region radius, forcing
% term and nodal weights), saved in RunInfo.Inverse.HessianVectorBasedState, provided it is for the same number of parameters and
% the same treatment of the bounds. Otherwise, or if CtrlVar.Inverse.HessianVectorBased.InitialTrustRegionRadius is given, the
% initial values are used. Previously each restart began with Delta=sqrt(number of parameters), which is typically much larger
% than the radius the run had settled at, so that the first steps of a restart run were rejected.
Restored="" ;
if isfield(RunInfo.Inverse,"HessianVectorBasedState") && isstruct(RunInfo.Inverse.HessianVectorBasedState)
    St=RunInfo.Inverse.HessianVectorBasedState ;
    if isfield(St,"Np") && St.Np==n && isfield(St,"Bounds") && St.Bounds==string(Settings.Bounds)
        if ~UserGivenRadius && isfinite(St.Delta) && St.Delta>0 ; Delta=St.Delta ; Restored=Restored+" Delta" ; end
        if isfinite(St.eta) && St.eta>0 ; eta=St.eta ; Restored=Restored+" eta" ; end
        if Settings.NodalTrustRegionScaling && numel(St.w)==n ; w=St.w(:) ; Restored=Restored+" nodal-weights" ; end
    end
end
if Restored~=""
    fprintf(" +++++++++++ Continuing with the trust-region state of the previous run (%s): Delta=%-g, eta=%-g\n",strtrim(Restored),Delta,eta)
else
    fprintf(" +++++++++++ Initial trust-region radius Delta=%-g\n",Delta)
end

for it=1:CtrlVar.Inverse.Iterations

    %% approximate solution of the Newton system, Steihaug-Toint
    Hv=HessianSetup(p,Fx,fOuts);

    if any(w>0)
        M=G+spdiags(w,0,n,n) ; dM=decomposition(M) ; PrecM=@(r) dM\r ;
        Mkey=[sum(w) sum(w.*(1:n)')] ;     % identifies M, for reusing the restricted factorisation
    else
        M=G ; PrecM=Prec ;
        Mkey=0 ;
    end

    if Projected
        [s,pred,cgInfo,Hs,StepCache]=ProjectedNewtonStep(p,g,Hv,M,PrecM,plb,pub,Delta,eta,Settings.MaxInnerIterations,StepCache,Mkey);
    else
        [s,Hs,cgInfo]=SteihaugToint(g,Hv,PrecM,M,Delta,eta*gM,Settings.MaxInnerIterations);
        pred=-(g'*s+0.5*(s'*Hs));     % predicted decrease
    end
    sG=sqrt(max(s'*(M*s),0));      % step length in the trust-region norm

    pt=p+s ;
    if Projected ; pt=min(max(pt,plb),pub) ; end   % (s is feasible by construction, this only removes round-off)
    [Jt,gt,~,fOutst,Fxt]=func(pt);
    ared=J-Jt;                    % actual decrease
    rho=ared/pred;
    if ~isfinite(Jt) ; rho=-Inf ; end

    %% trust-region update
    % mismatch between the actual and the predicted change of the gradient, and whether it is localised to a few variables
    Localised=false ;
    if Settings.NodalTrustRegionScaling && isfinite(Jt) && all(isfinite(gt))
        e2=(gt(:)-g-Hs).^2 ;
        [e2s,idx]=sort(e2,"descend") ;
        nTop=max(1,ceil(Settings.LocalisationFraction*n)) ;
        Localised=sum(e2s(1:nTop))>Settings.LocalisationShare*sum(e2s) ;
    end

    if rho<0.25 && Localised
        % the model failed, and the failure is localised: increase the weights of the variables concerned so that their step alone
        % would fill the trust region (at least four times the previous weight), but do not reduce the trust region elsewhere
        % The new weight is such that the step of each of these variables alone would fill the trust region (and at least four times
        % the previous weight), limited to 1e6 times the diagonal of G.
        I=idx(1:nTop) ;
        dG=full(diag(G)) ;
        sI=max(abs(s(I)),1e-3*max(abs(s))) ;
        w(I)=min(max(4*w(I),Delta^2./max(sI.^2,realmin)),1e6*dG(I)) ;
        Delta=0.5*Delta ;
    elseif rho<0.25
        Delta=0.25*sG;
    elseif rho>0.75 && (sG>=0.99*Delta || startsWith(cgInfo.Exit,"boundary") || startsWith(cgInfo.Exit,"neg. curv."))
        % the step was limited by the trust region (in the projected case the CG step in the free variables reaches the boundary
        % of the trust region remaining after the Cauchy step, while ||s||_G itself can be smaller than Delta)
        Delta=2*Delta;
    end
    if rho>0.75
        w=w/2 ;   % relax the nodal weights after good steps
        w(w<1e-12*max(full(diag(G))))=0 ;
    end

    Accepted=rho>1e-4;
    if Accepted
        gMold=gM ; [gModel,ModelNormCache]=GradientNorm(pt,g+Hs,G,Prec,plb,pub,Projected,ModelNormCache) ;   % norm of the gradient predicted by the model
        p=pt ; J=Jt ; g=gt(:) ; fOuts=fOutst ; Fx=Fxt ;
        [gM,NormCache]=GradientNorm(p,g,G,Prec,plb,pub,Projected,NormCache);
        etaNew=abs(gM-gModel)/max(gMold,realmin) ;                         % Eisenstat-Walker choice 1
        phi=(1+sqrt(5))/2 ;
        if eta^phi>0.1 ; etaNew=max(etaNew,eta^phi) ; end                  % safeguard
        eta=min(max(etaNew,1e-6),Settings.MaxForcingTerm) ;
        nRejected=0;
    else
        nRejected=nRejected+1;
    end

    [~,~,JBobs]=IndividualCostTerms(fOuts);
    fprintf("%5i  %5i  %-10s %15.10g  %15.10g  %15.10g  %11.4g  %11.4g  %11.4g  %9.4f",It0+it,cgInfo.Iterations,cgInfo.Exit,J,fOuts.I,fOuts.R,gM,sG,Delta,rho)
    if isBinversion ; fprintf("  %15.10g",JBobs) ; end
    if Projected ; fprintf("  %10i",sum(p<=plb | p>=pub)) ; end
    if Settings.NodalTrustRegionScaling ; fprintf("  %8i",sum(w>0)) ; end
    fprintf("  %8.1f s",toc(tStart)) ;   % cumulative wall-clock time
    if ~Accepted ; fprintf("  (step rejected)") ; end
    fprintf("\n")

    kNext=find(~isnan(RunInfo.Inverse.J),1,"last") ; if isempty(kNext) ; kNext=0 ; end
    RunInfo=StoreIteration(RunInfo,kNext+1,It0+it,J,fOuts,sG*Accepted,gM,cgInfo.Iterations,pred,Delta,Fields);
    % (16 Nov 2026) trust-region state, for continuing in a restart run (see above)
    RunInfo.Inverse.HessianVectorBasedState=struct("Delta",Delta,"eta",eta,"w",w,"Np",n,"Bounds",string(Settings.Bounds),"Iteration",It0+it) ;

    %% exit criteria
    if gM<=1e-10*gM0
        fprintf("UaOptimisationHessianVectorBased: gradient norm reduced by a factor of 1e10. Exiting.\n") ; break
    end
    if Delta<1e-10*Delta0 || nRejected>=20
        fprintf("UaOptimisationHessianVectorBased: trust region collapsed (no acceptable step found). Exiting.\n") ; break
    end

end

fprintf("\n +++++++++++ At end of inversion (HessianVectorBased):  J=%-g  I=%-g  R=%-g  |g|_G=%-g \n",J,fOuts.I,fOuts.R,gM)
if Projected
    fprintf(" +++++++++++ factorisations of G restricted to the free variables: %i (CG steps), %i (gradient norms)\n",StepCache.nFactorisations,NormCache.nFactorisations+ModelNormCache.nFactorisations)
end

end


function [s,Hs,Info]=SteihaugToint(g,Hv,Prec,G,Delta,tol,MaxIt)

% Approximate solution of min g's+1/2 s'Hs subject to ||s||_G <= Delta, preconditioned CG with the preconditioner G (Steihaug-Toint).
% Also returns H*s (accumulated, no additional Hessian-vector product needed).

s=zeros(size(g)) ; Hs=s ;
r=g ; z=Prec(r) ; d=-z ; rz=r'*z ;
Info.Iterations=0 ; Info.Exit="converged" ;

if sqrt(max(rz,0))<=tol ; return ; end

for k=1:MaxIt

    Info.Iterations=k ;
    Hd=Hv(d) ; dHd=d'*Hd ;
    Gs=G*s ; Gd=G*d ; sGs=s'*Gs ; sGd=s'*Gd ; dGd=d'*Gd ;

    if dHd<=0
        tau=ToBoundary(sGs,sGd,dGd,Delta) ;
        s=s+tau*d ; Hs=Hs+tau*Hd ; Info.Exit="neg. curv." ; return
    end

    alpha=rz/dHd ;
    if sGs+2*alpha*sGd+alpha^2*dGd>=Delta^2
        tau=ToBoundary(sGs,sGd,dGd,Delta) ;
        s=s+tau*d ; Hs=Hs+tau*Hd ; Info.Exit="boundary" ; return
    end

    s=s+alpha*d ; Hs=Hs+alpha*Hd ; r=r+alpha*Hd ;
    z=Prec(r) ; rzNew=r'*z ;

    if sqrt(max(rzNew,0))<=tol ; Info.Exit="converged" ; return ; end

    beta=rzNew/rz ; d=-z+beta*d ; rz=rzNew ;

end

Info.Exit="max. its" ;

end


function [s,pred,Info,Hs,Cache]=ProjectedNewtonStep(x,g,Hv,G,PrecG,l,u,Delta,eta,MaxInner,Cache,Gkey)

% Step for the bound-constrained trust-region subproblem, min q(s)=g's+1/2 s'Hs subject to l<=x+s<=u and ||s||_G<=Delta:
% (1) Cauchy step along the projected steepest-descent path, (2) Steihaug-Toint CG in the free variables, (3) projected search.

n=numel(x) ;
q=@(s,Hs) g'*s+0.5*(s'*Hs) ;

%% (1) Cauchy step: s(t)=P(x-t g)-x, backtracking from the minimiser of the model along -g (or the trust-region boundary)
Hg=Hv(g) ; gHg=g'*Hg ; gGg=g'*(G*g) ;
t=Delta/sqrt(max(gGg,realmin)) ;
if gHg>0 ; t=min(t,(g'*g)/gHg) ; end
for k=1:50
    sc=min(max(x-t*g,l),u)-x ;
    Hsc=Hv(sc) ;
    if sqrt(max(sc'*(G*sc),0))<=Delta*(1+1e-10) && q(sc,Hsc)<=0.01*(g'*sc)
        break
    end
    t=t/2 ;
end

xc=x+sc ; qc=q(sc,Hsc) ;
Active=(xc<=l) | (xc>=u) ;   % variables at their bounds at the Cauchy point
Free=~Active ;

s=sc ; Hs=Hsc ;
Info.Iterations=0 ; Info.Exit="Cauchy" ;

%% (2) Steihaug-Toint CG in the free variables, starting from the Cauchy point
DeltaW=Delta-sqrt(max(sc'*(G*sc),0)) ;     % so that ||sc+w||_G <= Delta
if any(Free) && DeltaW>0

    gc=g+Hsc ;                              % gradient of the model at the Cauchy point
    if all(Free) ; GF=G ; else ; GF=G(Free,Free) ; end
    [PrecF,Cache]=RestrictedSolver(Cache,G,PrecG,Free,Gkey) ;   % factorisation reused while Free (and G) do not change
    HvF=@(v) RestrictedHv(Hv,v,Free,n) ;
    gcF=gc(Free) ;
    gcM=sqrt(max(gcF'*PrecF(gcF),0)) ;

    if gcM>0

        [w,~,cgInfo]=SteihaugToint(gcF,HvF,PrecF,GF,DeltaW,eta*gcM,MaxInner) ;
        Info.Iterations=cgInfo.Iterations ; Info.Exit=cgInfo.Exit ;
        wFull=zeros(n,1) ; wFull(Free)=w ;

        %% (3) projected search along w: accept if the model decrease is at least that of the Cauchy step
        beta=1 ; Found=false ;
        for k=1:20
            st=min(max(xc+beta*wFull,l),u)-x ;
            Hst=Hv(st) ;
            if q(st,Hst)<=qc
                s=st ; Hs=Hst ; Found=true ;
                break
            end
            beta=beta/2 ;
        end
        if ~Found
            Info.Exit="Cauchy" ;
        elseif beta<1
            Info.Exit=Info.Exit+"*" ;  % projected search shortened the step
        end

    end

end

pred=-q(s,Hs) ;   % predicted decrease

end


function y=RestrictedHv(Hv,v,Free,n)

% Hessian-vector product restricted to the free variables

vFull=zeros(n,1) ; vFull(Free)=v ;
y=Hv(vFull) ;
y=y(Free) ;

end


function [gM,Cache]=GradientNorm(x,g,G,Prec,l,u,Projected,Cache)

% Gradient norm in the dual (inverse G) norm. For bound-constrained problems, the projected gradient: variables at a bound with the
% gradient pointing out of the feasible region are excluded.

if ~Projected
    gM=sqrt(max(g'*Prec(g),0)) ;
    return
end

Free=~((x<=l & g>0) | (x>=u & g<0)) ;
gF=g(Free) ;
if isempty(gF)
    gM=0 ;
else
    [SolveF,Cache]=RestrictedSolver(Cache,G,Prec,Free,0) ;   % factorisation reused while Free does not change (16 Nov 2026)
    gM=sqrt(max(gF'*SolveF(gF),0)) ;
end

end


function [Solve,Cache]=RestrictedSolver(Cache,G,PrecFull,Free,Key)

% (16 Nov 2026) Returns a function handle solving systems with G(Free,Free). If all variables are free, the given solver for G
% (PrecFull, eg the factorisation MUA.dG) is used. Otherwise the factorisation of G(Free,Free) is kept in Cache, and reused as long
% as Free and Key (which identifies G) are unchanged; previously G(Free,Free) was factorised in every iteration, and again for each
% (projected) gradient norm.

if all(Free)
    Solve=PrecFull ;
    return
end

if ~(isequal(Cache.Free,Free) && isequal(Cache.Key,Key) && ~isempty(Cache.dF))
    Cache.dF=decomposition(G(Free,Free)) ;
    Cache.Free=Free ;
    Cache.Key=Key ;
    Cache.nFactorisations=Cache.nFactorisations+1 ;
end
dF=Cache.dF ;
Solve=@(r) dF\r ;

end


function tau=ToBoundary(sGs,sGd,dGd,Delta)

% positive root of dGd tau^2 + 2 sGd tau + sGs - Delta^2 = 0, ie ||s+tau d||_G = Delta

tau=(-sGd+sqrt(max(sGd^2+dGd*(Delta^2-sGs),0)))/dGd ;

end


function RunInfo=StoreIteration(RunInfo,k,Iteration,J,fOuts,StepSize,GradNorm,InnerIterations,Decrement,TrustRegionRadius,Fields)

if numel(RunInfo.Inverse.J)<k
    Padding=nan(100,1);
    for ff=Fields ; RunInfo.Inverse.(ff)=[RunInfo.Inverse.(ff)(:) ; Padding] ; end
end

[Juv,Jhdot,JBobs]=IndividualCostTerms(fOuts);
RunInfo.Inverse.Iterations(k)=Iteration;
RunInfo.Inverse.J(k)=J;
RunInfo.Inverse.R(k)=fOuts.R;
RunInfo.Inverse.I(k)=fOuts.I;
RunInfo.Inverse.StepSize(k)=StepSize;
RunInfo.Inverse.GradNorm(k)=GradNorm;
RunInfo.Inverse.Decrement(k)=Decrement;
RunInfo.Inverse.Juv(k)=Juv;
RunInfo.Inverse.Jhdot(k)=Jhdot;
RunInfo.Inverse.JBobs(k)=JBobs;
RunInfo.Inverse.InnerIterations(k)=InnerIterations;
RunInfo.Inverse.TrustRegionRadius(k)=TrustRegionRadius;   % trust-region radius after the update in this iteration (16 Nov 2026)

end


function [Juv,Jhdot,JBobs]=IndividualCostTerms(fOuts)

Juv=NaN ; Jhdot=NaN ; JBobs=NaN ;
if isstruct(fOuts) && isfield(fOuts,"MisfitInfo") && isstruct(fOuts.MisfitInfo)
    if isfield(fOuts.MisfitInfo,"Juv") ; Juv=fOuts.MisfitInfo.Juv ; end
    if isfield(fOuts.MisfitInfo,"Jhdot") ; Jhdot=fOuts.MisfitInfo.Jhdot ; end
end
if isstruct(fOuts) && isfield(fOuts,"RegularisationInfo") && isstruct(fOuts.RegularisationInfo)
    if isfield(fOuts.RegularisationInfo,"JBobs") ; JBobs=fOuts.RegularisationInfo.JBobs ; end
end

end


function T=HVPType(CtrlVar)

T="-GN-";
if isfield(CtrlVar.Inverse,"HVP") && isfield(CtrlVar.Inverse.HVP,"Type") ; T=CtrlVar.Inverse.HVP.Type ; end

end
