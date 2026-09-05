






function [p,UserVar,RunInfo]=UaOptimisationHessianEstimate(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub)

narginchk(8,8)


%% Does an inversion using the Hessian
%
% Here the Hessian itself is constructed, i.e. not build iteratively as done in the BFSG method. 
%
% The Hessian is either approximated using brute-force finite-differences.
%
% Or by using the Direct-Adjoint method.
%
% This approach is memory hungry and building the Hessian can take a long time. 
%
% Consequently, this approach can only be used for small to medium sized problems with nodes on the order of few 10,000 or 
%
% The minimum is found using a two-dimensional (exact) trust-region approach.
%
% It is also possible to activate a Newton backtracking step, but this is mostly for testing purposes and once full confidence
% has been obtained in the trust-region approach, the Newton set can, and should, be deactivated. 
%
%
%%

%%
% good summary at
%
% https://uk.mathworks.com/help/releases/R2025b/optim/ug/constrained-nonlinear-optimization-algorithms.html#briahj8
%
% https://ecommons.cornell.edu/server/api/core/bitstreams/e68cc955-5b7c-4128-8707-522eb31ad378/content
%
% https://nmayorov.wordpress.com/2015/06/19/trust-region-reflective-algorithm/
%
% https://nmayorov.wordpress.com/2015/06/19/trust-region-reflective-algorithm/#more-207
%
% https://ecommons.cornell.edu/server/api/core/bitstreams/8830d99b-7b61-4271-ac76-f62504405051/content
%
% https://cs.uwaterloo.ca/~yuying/papers/Branch.pdf
%
% https://www.numerical.rl.ac.uk/media/people/nick-gould/ConnGoulToin88_mc.pdf
%
%
% https://math.stackexchange.com/questions/2814292/minimizing-convex-quadratic-with-box-constraints  (simple nice example of
% not just projecting onto the feasible domain)
%
% file:///C:/Users/Hilmar/Downloads/1984_Dembo_Tulowitzki_BQP_Yale_CS_tr302.pdf
%
%
% https://ccsp.hms.harvard.edu/wp-content/uploads/2022/10/Froehlich-Sorger-2022-Fides-Reliable-trust-region-optimization-for-parameter-estimation-of-ordinary-differential-equation-models.pdf
%
%%

%%
[isA,isB,isC] = isABC(CtrlVar);


if isA
    if CtrlVar.Inverse.Matern.logAGlen.alpha==1
     error("UaOptimisationHessianEstimate:WrongInputs","alpha=1 for the Matern parameter is not supported ")
    end
end


if isB
    if CtrlVar.Inverse.Matern.logB.alpha==1
        error("UaOptimisationHessianEstimate:WrongInputs","alpha=1 for the Matern parameter is not supported ")
    end
end


if isC
    if CtrlVar.Inverse.Matern.logC.alpha==1
        error("UaOptimisationHessianEstimate:WrongInputs","alpha=1 for the Matern parameter is not supported ")
    end
end


%%

doNewton=true;
doSteepestDescent=true; 
doTrustRegion=true;


% Tolerances: (me to do)  These are currently hardwired, but should be included in CtrlVar.
SubOptimalityTolerance=1e-10;
JTolerance=0.01;
dJTolerance=0.0 ; 
dpTolerance=0.0;


%%
p=p0;


MaxNewtonSteps=CtrlVar.Inverse.Iterations;

Jvector=nan(MaxNewtonSteps+1,1);
SlopeVector=nan(MaxNewtonSteps+1,1);
gammaVector=nan(MaxNewtonSteps+1,1);
GradNormVector=nan(MaxNewtonSteps+1,1);
SubOptimalityVector=nan(MaxNewtonSteps+1,1);
nPar=numel(p);

iRange=1:nPar; % all of them

Delta_new=nan;
iIteration=0; lStart=0 ; 
gammaSteepestDescentLast=nan;
JTrustRegion=inf ; DeltaMax=nan;
TrustRegionStepAccepted=false;


%% Build metric matrix, i.e. the Gram matrix


if  ~isfield(MUA,"MetricMatrix") || isempty(MUA.MetricMatrix)
    error("UaOptimisationHessianEstimate:MetricMatrix","Was expecting MetricMatrix to be a field of MUA, but it is not. \n")
end


MetricMatrix=MUA.MetricMatrix;
dMetricMatrix=decomposition(MetricMatrix);


while true

    iIteration=iIteration+1;

    if contains(CtrlVar.Inverse.MinimisationMethod,"BruteForceHessian")

        [~,Hessian,g0,J0] = CalcBruteForceHessian(func,p,CtrlVar,iRange) ;

        if isnan(J0)
            error("UaOptimisationHessianEstimate:J0IsNaN","NaN in J0")
        end

    elseif contains(CtrlVar.Inverse.MinimisationMethod,"DirectAdjointHessian")

        [J0,g0,Hessian]=func(p);

        if CtrlVar.Inverse.TestDirectAdjoint.isTrue

            %%
            iCol=randi(numel(p));
            % iCol=88;
            [~,HessianFD,g0,J0] = CalcBruteForceHessian(func,p,CtrlVar,iCol) ;
            Diff=norm(HessianFD(:,iCol)-Hessian(:,iCol))/norm(Hessian(:,iCol));
            fprintf("UaOptimisationHessianEstimate: normalised norm of difference between Direct-Adjoint and finite-difference Hessian for column %i is: %g \n",iCol,Diff)
            FigHessDA_FD=FindOrCreateFigure("Test: DirectAdjoint Hess") ; clf(FigHessDA_FD)
            hold off ;
            plot(HessianFD(:,iCol),Hessian(:,iCol),"o") ;
            axis equal ;
            AX=axis;
            hold on ;
            plot([AX(1) AX(2)],[AX(1) AX(2)],"--")
            title(sprintf("FD test of Direct-Adjoint Hessian for column %i",iCol),Interpreter="latex",FontWeight="bold",FontSize=14)
            subtitle(sprintf("Inverting for: %s ",CtrlVar.Inverse.InvertFor),Interpreter="latex")
            xlabel("Finite Differences",Interpreter="latex")
            ylabel("Direct Adjoint",Interpreter="latex")
            ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off
            %%


        end

    end



    [Hessian,lStart]=CheckIfHessianIsSPDandIfNotMakeItSo(Hessian,MetricMatrix,lStart) ;
    lCondition=1e-5; lConditionMin=0;
    Hessian=ImproveMatrixCondition(Hessian,MetricMatrix,lCondition,lConditionMin) ;

    dpNewton=Hessian\(-g0);  % Here I need to add in the BCs, I need BCs on dp, i.e. dA and dC
    % To do: Currently there are not BCs applied to the inversion fields, but this should be an option going forward.
    slope0Newton=g0'*dpNewton;

    if anynan(dpNewton)
        fprintf("Solving the Newton system resulted in nan. \n")
        error("UaOptimisationHessianEstimate:dpIsNaN","NaN in dp")
    end

    if anynan(p)
        fprintf("p contains nan. \n")
        error("UaOptimisationHessianEstimate:pIsNaN","NaN in p")
    end

    if slope0Newton >0
        fprintf("Slope at origin in Newton direction is positive! \n")
        gammaNewton=nan;
        JNewton=inf;
  
    end

    if doNewton && slope0Newton < 0

        gamma=1;  % This is the initial guess for a good gamma, now that the DA Hessian is accurate, seems best to use this by default
        J1=func(p+gamma*dpNewton);
        while isnan(J1)
            gamma=gamma/10;
            J1=func(p+gamma*dpNewton);
        end
        Func=@(gamma) func(p+gamma*dpNewton); % here a plus sign because I'm going in the direction dp, this is the Newton direction

        CtrlVar.NewtonAcceptRatio=0.1 ;CtrlVar.BacktrackingGammaMin=gamma/1e5; CtrlVar.LineSearchAllowedToUseExtrapolation=false;
        CtrlVar.InfoLevelBackTrack=100 ; CtrlVar.doplots=1 ;
        [gammaNewton,JNewton]=BackTracking(slope0Newton,gamma,J0,J1,Func,CtrlVar);

    end




    % To do: Here I must add the boundary conditions for A/B/C. At the moment I do not prescribe any BCs for any of these fields
    % so there is nothing to add here currently. However, I intent to add in the boundary conditions for the B inversion, but the
    % B inversion is currently under development.
    dpSteepestDescent=dMetricMatrix\(-g0); % pre-multiplying, note that I must use the inverse...!
    slope0SteepestDescent=g0'*dpSteepestDescent;

    if slope0SteepestDescent >0
        fprintf("Slope at origin in steepest descent direction is positive! \n")
        gammaSteepestDescent=nan;
        JSteepestDescent=inf;
        doSteepestDescent=0;
    end

    if doSteepestDescent

        % If I believe in the quadratic model then
        %
        % m(p)=f + g p + 0.5 p H p
        %
        % line search: g p \gamma +0.5 p H p \gamma^2
        %
        % If I select the search direction p as p= -G\g
        % 
        % minimum found for:
        %
        % $$ \gamma=-\frac{g^T p}{p^T H p} $$
        % 
        % If $G=H$ then $p=-H \ g$ that is $g=-H p$ 
        % I have the newton step and the best gamma is $\gamma=1$
        %
        %
        %
        %

        %if isnan(gammaSteepestDescentLast)
        % Estimate a good gamma from the minimum along the search direction
        gamma=- g0' * dpSteepestDescent/(dpSteepestDescent'*Hessian*dpSteepestDescent);
        % else
        %     % Not sure I always trust the quadratic model that well. Maybe best to estimate a reasonable guess for gamma based on last
        %     % gamma min value, increased by a factor
        %     gamma=2*gammaSteepestDescentLast;
        % end

        J1=func(p+gamma*dpSteepestDescent);
        while isnan(J1)
            gamma=gamma/10;
            J1=func(p+gamma*dpSteepestDescent);
        end
        Func=@(gamma) func(p+gamma*dpSteepestDescent); % here a plus sign because I'm going in the direction dp, this is the Newton direction

        CtrlVar.NewtonAcceptRatio=0.1 ;CtrlVar.BacktrackingGammaMin=gamma/1e5; CtrlVar.LineSearchAllowedToUseExtrapolation=false;
        CtrlVar.InfoLevelBackTrack=100 ; CtrlVar.doplots=1 ;
        [gammaSteepestDescent,JSteepestDescent]=BackTracking(slope0SteepestDescent,gamma,J0,J1,Func,CtrlVar);
        gammaSteepestDescentLast=gammaSteepestDescent;

    end


   

    %% testing two-dimensional subspace approach
    if doTrustRegion

        if isnan(DeltaMax)

            [isA,isB,isC] = isABC(CtrlVar);

            Neff=0;

            if isA
                Neff  = Neff+ MUA.Area/CtrlVar.Inverse.Matern.logAGlen.rho^2;
            end

            if isB
                Neff  = Neff + MUA.Area/CtrlVar.Inverse.Matern.B.rho^2;
            end

            if isC
                Neff  = Neff + MUA.Area/CtrlVar.Inverse.Matern.logC.rho^2;
            end




            DeltaRef = 0.8* CtrlVar.TrustRegion.nSigma*sqrt(Neff);

            % DeltaMax = min(10*sqrt(dpNewton'*MetricMatrix*dpNewton), 10*DeltaRef);
            % DeltaMax = max(DeltaMax, DeltaRef/10);

            DeltaMax=DeltaRef;
            DeltaMin=DeltaMax/1e6;

            % The way DeltaMax is selected, only works for alpha>1 in the Matern formulation. This will also not work if one is using the old Tikhonov approach.
            % going forward, I'm simply going to insist on using alpha>1.
            %
            % The other inversion parts of the code that do not use the subspace minimization still work for alpha=1.
            %
            assert(isfinite(DeltaMax) && DeltaMax>0, ...
                "UaOptimisationHessianEstimate:BadDeltaMax","DeltaMax=%g",DeltaMax)

        end


        if isnan(Delta_new)
            Delta=DeltaMax;
        else
            Delta=Delta_new;
        end

        TrustRegionStepAccepted=false;

        while ~TrustRegionStepAccepted && Delta>DeltaMin

            [dpTR,yTR,info]=TwoDSubspaceTrustRegionGram(g0,Hessian,dpNewton,dpSteepestDescent,Delta,MetricMatrix) ;
            JTrustRegion=func(p+dpTR)   ;
            ActualReduction = -(JTrustRegion-J0) ; % actual reduction for the
            PredicedReduction = -(g0'*dpTR + 0.5 * dpTR' * Hessian * dpTR) ; % predicted reduction (note that J0 cancels)
            rho=ActualReduction/PredicedReduction ;

            [p0_new, Delta_new, TrustRegionStepAccepted] = TrustRegionUpdate(p, dpTR, J0, JTrustRegion, PredicedReduction, Delta, DeltaMax,info.stepnormG) ;

            fprintf("TrustRegion: accepted=%s \t case=%s \t rho=%f \t Delta=%f \t delta_new=%g  dp=%3.3f dNewton  %+3.3f dSteepestDescent \n ",string(TrustRegionStepAccepted),info.case,rho,Delta,Delta_new,yTR(1),yTR(2))



            Delta=Delta_new;

        end

        if Delta==DeltaMax && rho >0.95 && rho < 1.05
            fprintf("Delta is hitting against DeltaMax for rho=%g\n",rho)
            fprintf("Note: If Delta sits at DeltaMax repeatedly while steps are accepted with rho near 1, the cap is binding and should be raised. \n")
            fprintf("      The place to change this is by increasing the value of CtrlVar.TrustRegion.nSigma, currently at %f\n",CtrlVar.TrustRegion.nSigma)
        end

    end

    %%
    % gammaSteepestDescent=nan; JSteepestDescent=nan ; gammaSDMax=nan 
    
    gammaNewtonMax=1.5; gammaSDMax=gammaSteepestDescent*2;
    PlotCostVersusStepSizeAlongNewtonDirection(func,p,dpNewton,g0,Hessian,gammaNewton,JNewton,dpSteepestDescent,gammaSteepestDescent,JSteepestDescent,gammaNewtonMax,gammaSDMax,doSteepestDescent);
    drawnow

    fprintf("====> JNewton/J0=%g \t JSteepestDescent/J0=%g \t JTrustRegion/J0=%g \n",JNewton/J0,JSteepestDescent/J0,JTrustRegion/J0)

    if doNewton
        rhoJNewton=JNewton/J0;
    else
        rhoJNewton=inf;
    end
 

    if doSteepestDescent
        rhoJSteepestDescent=JSteepestDescent/J0;
    else
        rhoJSteepestDescent=inf;
    end


    if doTrustRegion && TrustRegionStepAccepted
        rhoJTrustRegion=JTrustRegion/J0;
    else
        rhoJTrustRegion=inf;
    end

    % pick the best
    [rhoJmin,iJmin]=min([rhoJNewton,rhoJSteepestDescent,rhoJTrustRegion]);

    if rhoJmin> 1

        fprintf(" Value of cost function could not be reduced. \n")
        break

    end

    switch iJmin

        case 1
            Direction="Newton";
        case 2
            Direction="SteepestDescent";
        case 3
            Direction="TrustRegion";
    end

    switch Direction

        case "Newton"


            dp=gammaNewton*dpNewton;
            p=p+dp;
         
            dpNorm=norm(dp)/norm(p);
            J=JNewton;
            gamma=gammaNewton;

        case "SteepestDescent"


            dp=gammaSteepestDescent*dpSteepestDescent;
            p=p+dp;
          
            dpNorm=norm(dp)/norm(p);
            J=JSteepestDescent;
            gamma=gammaSteepestDescent;

        case "TrustRegion"

            p=p0_new;

            dpNorm=norm(dpTR)/norm(p);
            J=JTrustRegion;
            gamma=nan;
    end

   SubOptimality=-g0'*dpNewton/2  ; % Newton decrement g0' H^{-1} g /2
   

    dJ=J0-J;


    pUpperViolation=pub-p;
    if any(pUpperViolation<0)
        fprintf("p violates pub! \n")
    end

    fprintf("%3i:(%s) \t gamma=%g \t J=%g \t J0=%g \t sub-obtimality=%g \t |dp|/|p|=%g \t J/J0=%g \n",iIteration,Direction,gamma,J,J0,SubOptimality,dpNorm,J/J0)

    Jvector(iIteration)=J0;
    Jvector(iIteration+1)=J;
    SlopeVector(iIteration)=slope0Newton;
    gammaVector(iIteration)=gamma;
    SubOptimalityVector(iIteration)=SubOptimality;
    GradNormVector(iIteration)=norm(g0)/sqrt(numel(g0));


    drawnow

    % exit?
    if iIteration>(MaxNewtonSteps-1)
        fprintf("maximum number of iterations reached. \n")
        break
    end

    if slope0Newton< 0 && SubOptimality<SubOptimalityTolerance
        fprintf("subtolerance (%g) reached with %g. \n",SubOptimalityTolerance,SubOptimality)
        break
    end

    if dJ<dJTolerance
        fprintf("dJ tolerance (%g) reached with dJ=%g. \n",dJTolerance,dJ)
        break
    end

    if J<JTolerance
        fprintf("J tolerance (%g) reached with J=%g. \n",JTolerance,J)
        break
    end

    if dpNorm< dpTolerance
        fprintf("dp tolerance (%g) reached with dpNorm=%g. \n",dpTolerance,dpNorm)
        break
    end

    % if J/J0>=0.999999
    %     fprintf("stagnated. \n")
    %     break
    % end

end

I=~isnan(Jvector);
Jvector=Jvector(I);
SubOptimalityVector=SubOptimalityVector(I);
itVector=0:(numel(Jvector)-1) ; itVector=itVector(:);

I=~isnan(GradNormVector); GradNormVector=GradNormVector(I); GradNormVector=GradNormVector(:) ;  GradNormVector=[GradNormVector;NaN]; % Make sure it has the same length as itVector
I=~isnan(gammaVector); gammaVector=gammaVector(I); gammaVector=gammaVector(:) ;  gammaVector=[gammaVector;NaN]; % Make sure it has the same length as itVector

figIt=FindOrCreateFigure("J iteration") ; clf(figIt)
yyaxis left
semilogy(itVector,Jvector,"ob-",LineWidth=2,DisplayName="$J$",MarkerFaceColor="b") ;
hold on
semilogy(itVector,SubOptimalityVector,"sg-",DisplayName="sub-optimality")
ylabel("$J$",Interpreter="latex") ;
yyaxis right
hold off
semilogy(itVector,GradNormVector,"or-",LineWidth=2,DisplayName="$\|\nabla J \|$") ;
ylabel("$\|\nabla J \|$",Interpreter="latex") ;
xlabel("Iteration")
legend(Interpreter="latex");


itRestart=max(RunInfo.Inverse.Iterations);

RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;itVector+itRestart];
RunInfo.Inverse.J=[RunInfo.Inverse.J;Jvector];
RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;GradNormVector];
RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gammaVector];




end

%%