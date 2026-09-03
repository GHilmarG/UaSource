






function [p,UserVar,RunInfo]=UaOptimisationHessianEstimate(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub)

narginchk(8,8)

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



p=p0;


MaxNewtonSteps=CtrlVar.Inverse.Iterations;

Jvector=nan(MaxNewtonSteps+1,1);
SlopeVector=nan(MaxNewtonSteps+1,1);
gammaVector=nan(MaxNewtonSteps+1,1);
GradNormVector=nan(MaxNewtonSteps+1,1);
SubOptimalityVector=nan(MaxNewtonSteps+1,1);
nPar=numel(p);

iRange=1:nPar; % all of them




SubOptimalityTolerance=1e-10;
JTolerance=0.01;
dJTolerance=0.0;
dpTolerance=0.0;
Delta_new=nan;
iIteration=0; lStart=0 ; gammaSDLast=inf;

doNewton=true;
doSteepestDescent=false;
doTrustRegion=true;


JTrustRegion=inf ; DeltaMax=nan;
JNewton=inf;
JSteepestDescent=inf;

dpSteepestDescent=nan;
gammaSteepestDescent=nan;
gammaSDMax=nan;
gammaNewton=nan;
MetricMatrix=[];

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

    %% Build metric matrix, i.e. the Gram matrix
    if isempty(MetricMatrix)  % Only build it once, and reuse the decomposition

        MetricMatrix=BuildMetricMatrix(CtrlVar,MUA); 
        dMetricMatrix=decomposition(MetricMatrix);
    end
    %%
 

    [Hessian,lStart]=CheckIfHessianIsSPDandIfNotMakeItSo(Hessian,MetricMatrix,lStart) ;
    lCondition=1e-5; lConditionMin=0;
    Hessian=ImproveMatrixCondition(Hessian,MetricMatrix,lCondition,lConditionMin) ;

    dpNewton=Hessian\(-g0);  % Here I need to add in the BCs, I need BCs on dp, i.e. dA and dC
    slope0Newton=g0'*dpNewton;

    if anynan(dpNewton)
        fprintf("Solving the Newton system resulted in nan. \n")
        error("UaOptimisationHessianEstimate:dpIsNaN","NaN in dp")
    end

    if anynan(p)
        fprintf("p contains nan. \n")
        error("UaOptimisationHessianEstimate:pIsNaN","NaN in p")
    end


    CtrlVar.GradientReflective=false;

    if CtrlVar.GradientReflective

        % This just about OK.
        %
        % The Reflective Transformation, R(p), only reflects points that are outside of the box constraints. If for some i, p_i=l_i
        % or p_i=u_i, the value of p_i is not changed. This means that for a given search direction (-g) direction, g, the line-search may result
        % in values being exactly at the boundary, but never outside the box constraints.
        %
        % Once the gradient is recalculated, i.e. at the start of next iteration, the corresponding elements of the search direction
        % are flipped where p is at the boundary and the unmodified search direction points out. This ensures that the next update
        % will shift p away from the box boundary. However, now the gradient is no longer the gradient of the cost function!
        %
        % 1) change the sign of a gradient element whenever a box constraint is violated
        % 2)
        %
        %
        gammaUpperVector=(pub-p)./dpNewton;
        gammaUpperVector(gammaUpperVector<eps)=nan ;  % where this is negative, there is no constraint on the gamma
        [gammaNewtonMax,Imin]=min(gammaUpperVector)  ; % this is the smallest positive gamma that does not violate

        % What to do where p_i=pub_i and dp_i > 0  ? Then any finite positive step size will violate pub at those locations
        % This would cause zero step size, or more generally, a very small step size if p_i is very close to pub_i and dp_i >0.
        %
        % Here one can try 'reflection' where p_i is reflected by setting it to -p_i


        %% reflection
        gammaNewtonMin=0.2;
        if gammaNewtonMax < gammaNewtonMin

            I=find(gammaVector<gammaNewtonMin);
            dpNewton(I)=-dpNewton(I) ;      % Reflection

            gammaVector=(pub-p)./dpNewton;
            gammaVector(gammaVector<eps)=nan ;
            [gammaNewtonMax,Imin]=min(gammaVector)  ;

        end
    else
        gammaNewtonMax=inf;
    end
    %%

    if doNewton

        gamma=1;  % This is the initial guess for a good gamma, now that the DA Hessian is accurate, seems best to use this by default

        if gammaNewtonMax< gamma
            gamma=gammaNewtonMax;
            CtrlVar.LineSearchAllowedToUseExtrapolation=false;
            pUpperViolation=pub-(p+gamma*dpNewton) ;
            UaPlots(CtrlVar,MUA,[],pUpperViolation,FigureTitle="pUpperViolation Newton")  ;
            CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
            min(pUpperViolation)

        end



        J1=func(p+gamma*dpNewton);

        while isnan(J1)
            gamma=gamma/10;
            J1=func(p+gamma*dpNewton);
        end

        Func=@(gamma) func(p+gamma*dpNewton); % here a plus sign because I'm going in the direction dp, this is the Newton direction

        CtrlVar.InfoLevelBackTrack=0;
        % Newton direction
        slope0Newton=g0'*dpNewton;

        if slope0Newton<0
            CtrlVar.NewtonAcceptRatio=0.1 ;CtrlVar.BacktrackingGammaMin=gamma/1e5; CtrlVar.LineSearchAllowedToUseExtrapolation=false;
            CtrlVar.InfoLevelBackTrack=100 ; CtrlVar.doplots=1 ;
            [gammaNewton,JNewton]=BackTracking(slope0Newton,gamma,J0,J1,Func,CtrlVar);
        else
            fprintf("Slope at origin in Newton direction is positive! \n")
            gammaNewton=nan;
            JNewton=inf;
        end

    end


   


    dpSteepestDescent=dMetricMatrix\(-g0); % pre-multiplying, note that I must use the inverse...!

    % steepest decent
    if doSteepestDescent

        Func=@(gamma) func(p+gamma*dpSteepestDescent);

        slope0=g0'*dpSteepestDescent;

        if slope0 > 0
            fprintf("Slope at origin in steepest-descent direction is positive! \n")
        end

        % calculate maximum step-size that does not violate upper limit
        %if CtrlVar.GradientReflective

        gammaUpperVector=(pub-p)./dpSteepestDescent;
        gammaUpperVector(gammaUpperVector<eps)=nan ;  % where this is negative, there is no constraint on the gamma

        gammaSDMax=min(gammaUpperVector);



        gamma=min(2*gammaSDLast,gammaSDMax);
        J1=Func(gamma);
        while isnan(J1)
            gamma=gamma/10;
            J1=func(p+gamma*dpSteepestDescent);
        end

        gammaSDLast=min(gammaSDLast,gamma);

        CtrlVar.NewtonAcceptRatio=0.1 ;CtrlVar.BacktrackingGammaMin=gammaSDLast/1e6;
        [gammaSteepestDescent,JSteepestDescent]=BackTracking(slope0,gamma,J0,J1,Func,CtrlVar);
        gammaSDLast=gammaSteepestDescent;



    end



    %% testing two-dimensional subspace approach
    if doTrustRegion

        if isnan(DeltaMax)
            % only set DeltaMax once, and set it to a generously large value
            DeltaMax=10*sqrt(dpNewton'*MetricMatrix*dpNewton);
            DeltaMin=DeltaMax/1e6;
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

            [p0_new, Delta_new, TrustRegionStepAccepted] = TrustRegionUpdate(p, dpTR, J0, JTrustRegion, PredicedReduction, Delta, DeltaMax) ;

            fprintf("TrustRegion: accepted=%s \t case=%s \t rho=%f \t Delta=%f \t delta_new=%g  dp=%3.3f dNewton  %+3.3f dSteepestDescent \n ",string(TrustRegionStepAccepted),info.case,rho,Delta,Delta_new,yTR(1),yTR(2))

            % if Delta==DeltaMax
            %        fprintf("Delta is hitting against DeltaMax for rho=%g\n",rho)
            % end

            Delta=Delta_new;

        end

    end

    %%
    % PlotCostVersusStepSizeAlongNewtonDirection(func,p,dpNewton,g0,Hessian,gammaNewton,JNewton,dpSteepestDescent,gammaSteepestDescent,JSteepestDescent,gammaNewtonMax,gammaSDMax,doSteepestDescent);

    %fprintf("====> JNewton/J0=%g \t JSteepestDescent/J0=%g \t JTrustRegion/J0=%g \n",JNewton/J0,JSteepestDescent/J0,JTrustRegion/J0)

    rhoJNewton=JNewton/J0;
    rhoJSteepestDescent=JSteepestDescent/J0;
 
    if TrustRegionStepAccepted
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

    if SubOptimality<SubOptimalityTolerance
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

    if J/J0>=0.999999
        fprintf("stagnated. \n")
        break
    end

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
lg=legend(Interpreter="latex");


itRestart=max(RunInfo.Inverse.Iterations);

RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;itVector+itRestart];
RunInfo.Inverse.J=[RunInfo.Inverse.J;Jvector];
RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;GradNormVector];
RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gammaVector];




end

%%