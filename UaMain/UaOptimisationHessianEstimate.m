

function [p,UserVar,RunInfo]=UaOptimisationHessianEstimate(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub)

narginchk(8,8)

%%
%load("TestSaveH.mat","func","p0","CtrlVar","iRange","MUA","F")  ;


%%

p=p0;


MaxNewtonSteps=CtrlVar.Inverse.Iterations;

Jvector=nan(MaxNewtonSteps+1,1);
SlopeVector=nan(MaxNewtonSteps+1,1);
gammaVector=nan(MaxNewtonSteps+1,1);
GradNormVector=nan(MaxNewtonSteps+1,1);

nPar=numel(p);

iRange=1:nPar; % all of them




SubOptimalityTolerance=0;
dJTolerance=0.0;
dpTolerance=0.0;

iNewton=0; lStart=0 ; gammaSDLast=inf; gammaNewtonLast=1; 

while true

    iNewton=iNewton+1;

    if contains(CtrlVar.Inverse.MinimisationMethod,"BruteForceHessian")

        [Hsparse,Hfull,g0,J0] = CalcBruteForceHessian(func,p,CtrlVar,iRange) ;

        if isnan(J0)
            error("UaOptimisationHessianEstimate:J0IsNaN","NaN in J0")
        end

    elseif contains(CtrlVar.Inverse.MinimisationMethod,"DirectAdjointHessian")

        [J0,g0,Hfull]=func(p);

        if CtrlVar.Inverse.TestDirectAdjoint.isTrue

            %%
            iCol=randi(numel(p));
           % iCol=88;
            [~,HfullFD,g0,J0] = CalcBruteForceHessian(func,p,CtrlVar,iCol) ;  
            Diff=norm(HfullFD(:,iCol)-Hfull(:,iCol))/norm(Hfull(:,iCol));
            fprintf("UaOptimisationHessianEstimate: normalised norm of idfference between Direct-Adjoint and finite-difference Hessian for column %i is: %g \n",iCol,Diff)
            FigHessDA_FD=FindOrCreateFigure("Test: DirectAdjoint Hess") ; clf(FigHessDA_FD)
            hold off ; 
            plot(HfullFD(:,iCol),Hfull(:,iCol),"o") ; 
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




    [Hfull,lStart]=CheckIfHessianIsSPDandIfNotMakeItSo(Hfull,MUA,lStart) ;
    lCondition=1e-5; lConditionMin=0;
    Hfull=ImproveMatrixCondition(Hfull,MUA.M,lCondition,lConditionMin) ;

    dp=Hfull\(-g0);  % Here I need to add in the BCs, I need BCs on dp, i.e. dA and dC




    %   [L,cuv]=AssembleLuvSSTREAM(CtrlVar,MUA,BCs) ;
    %
    % L=chol(Hfull) ;
    % tol=1e-6; maxit=30 ;
    % [dpTest,fl,rr,it,rv1,rvcgl]=minres(Hfull,-g0,tol,maxit,L,L');
    % Fig=figure ; plot(0:length(rv1)-1,rv1/norm(g0),"-or") ; ax=gca ; ax.YScale="log";
    %
    % L=chol(Hfull) ;
    % tol=1e-6; maxit=30 ;
    %
    % afun=@(x) HVP(x,Hfull) ;
    % [dpTest,fl,rr,it,rv1,rvcgl]=minres(afun,-g0,tol,maxit,L,L');
    % FigHVP=figure ; plot(0:length(rv1)-1,rv1/norm(g0),"-or") ; ax=gca ; ax.YScale="log";
    %
    %
    %
    % function y=HVP(x,Hfull)
    %     y=Hfull*x;
    %
    % end
    % HVP=HessianVectorProduct(p,d,func)

    %    UaPlots(CtrlVar,MUA,[],dp,FigureTitle="Newton dp") ; CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);



    %  D=norm(dp) ;   H=Hfull ; l=0 ; g=g0;  E=blkdiag(MUA.M,MUA.M) ; l=TrustRegionSubproblem(H,E,g,l,D) ;


    if anynan(dp)
        fprintf("Solving the Newton system resulted in nan. \n")
        error("UaOptimisationHessianEstimate:dpIsNaN","NaN in dp")
    end

    if anynan(p)
        fprintf("p contains nan. \n")
        error("UaOptimisationHessianEstimate:pIsNaN","NaN in p")
    end




    % Since the Hessian may have been modified, it is not clear what a sensible first step could be.
    %
    % Also, the Hessian is not exact.
    %
    % Furthermore, what is the largest gamma I can use without violating the limits?
    %

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
        gammaUpperVector=(pub-p)./dp;
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
            dp(I)=-dp(I) ;      % Reflection

            gammaVector=(pub-p)./dp;
            gammaVector(gammaVector<eps)=nan ;
            [gammaNewtonMax,Imin]=min(gammaVector)  ;

        end
    else
        gammaNewtonMax=inf;
    end
    %%

    gamma=min(2*gammaNewtonLast,1);  %take last estimate for the min found in backtracking, and expand the radius by a factor of two each time.
    if gamma == 0 % this could have happened if previous line-search was not successful 
        gamma=1;
    end
    if gammaNewtonMax< gamma 
        gamma=gammaNewtonMax;
        CtrlVar.LineSearchAllowedToUseExtrapolation=false;
        pUpperViolation=pub-(p+gamma*dp) ;
        UaPlots(CtrlVar,MUA,[],pUpperViolation,FigureTitle="pUpperViolation Newton")  ;
        CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
        min(pUpperViolation)

    end


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

    % gamma=1;
    J1=func(p+gamma*dp);

    while isnan(J1)
        gamma=gamma/10;
        J1=func(p+gamma*dp);
    end

    Func=@(gamma) func(p+gamma*dp); % here a plus sign because I'm going in the direction dp, this is the Newton direction

    CtrlVar.InfoLevelBackTrack=0;
    % Newton direction
    slope0=g0'*dp;

    if slope0<0
        CtrlVar.NewtonAcceptRatio=0.1 ;CtrlVar.BacktrackingGammaMin=gamma/1e5; CtrlVar.LineSearchAllowedToUseExtrapolation=false;
         CtrlVar.InfoLevelBackTrack=100 ; CtrlVar.doplots=1 ;
        [gammaNewton,JNewton,BackTrackInfo]=BackTracking(slope0,gamma,J0,J1,Func,CtrlVar);
    else
        fprintf("Slope at origin in Newton direction is positive! \n")
        gammaNewton=nan;
        JNewton=inf;
    end

    gammaNewtonLast=gammaNewton;

    doSteepestDecent=true;
    if doSteepestDecent

        % I also do Steepest decent, and then compare.
        % This is because calculating the Hessian is so expensive and, in comparison, the line-search cheap

        % Then maybe modify H and solve until slope negative, or just go for the gradient direction
        nNodes=size(MUA.M,1);
        nH=size(Hfull,1);
        M=MUA.M ;
        D=MUA.Dxx+MUA.Dyy;
        
        if nH==nNodes
          
        elseif nH==2*nNodes
            M=blkdiag(M,M) ;
            D=blkdiag(D,D) ;

        else
            error("wrong dimensions")
        end

    %   CtrlVar.Inverse.AdjointGradientPreMultiplier="H1";


        if CtrlVar.Inverse.AdjointGradientPreMultiplier=="L2"
            
            P=M;

        elseif CtrlVar.Inverse.AdjointGradientPreMultiplier=="H1"

            ga=CtrlVar.Inverse.PreMultiplier.H1.ga;
            gs=CtrlVar.Inverse.PreMultiplier.H1.gs;
            P=gs*D+ga*M;

        else

            P=1;
        end

        g0SD=P\(-g0); % pre-multiplying, note that I must use the inverse...!


        Func=@(gamma) func(p+gamma*g0SD);

        slope0=g0'*g0SD;
        gammaSlope0=-0.1*J0/slope0;

        % calculate maximum step-size that does not violate upper limit
        %if CtrlVar.GradientReflective

        gammaUpperVector=(pub-p)./g0SD;
        gammaUpperVector(gammaUpperVector<eps)=nan ;  % where this is negative, there is no constraint on the gamma

        gammaSDMax=min(gammaUpperVector);

        %else
        %    gammaSDMax=inf;
        %end

        % gamma=max(gammaSlope0,gammaSDMax);

        % gamma=gammaSDMax ; pUpperViolation=pub-(p+gamma*g0SD) ;  UaPlots(CtrlVar,MUA,[],pUpperViolation,FigureTitle="pUpperViolation SD")  ;  CM=cmocean('balanced',25,'pivot',0) ; colormap(CM); min(pUpperViolation)



        % J1=Func(gamma);
        
        gamma=min(2*gammaSDLast,gammaSDMax);   
        J1=Func(gamma);
        while isnan(J1)
            gamma=gamma/10;
            J1=func(p+gamma*g0SD);
        end

        gammaSDLast=min(gammaSDLast,gamma);

        CtrlVar.NewtonAcceptRatio=0.1 ;CtrlVar.BacktrackingGammaMin=gammaSDLast/1e6;
        [gammaSD,JSD]=BackTracking(slope0,gamma,J0,J1,Func,CtrlVar);
        gammaSDLast=gammaSD;

    else

        JSD=inf;
        gammaSD=nan;
        g0SD=nan; 
        gammaSDMax=nan;

    end

    fprintf("====> JNewton/J0=%g \t JSD/J=%g \n",JNewton/J0,JSD/J0)

    if JNewton  < JSD
        fprintf("Newton step resulted in a greater reduction than steepest decent \n ")
        SubOptimality=-g0'*dp/2  ;
        dpNorm=norm(gamma*dp)/norm(p);
        J=JNewton;
        gamma=gammaNewton;
        Direction="Newton";
    else
        SubOptimality=inf;
        dpNorm=inf;
        J=JSD;
        gamma=gammaSD;
        Direction="gradient";
    end



    PlotCostVersusStepSizeAlongNewtonDirection(func,p,dp,g0,gammaNewton,JNewton,g0SD,gammaSD,JSD,gammaNewtonMax,gammaSDMax,doSteepestDecent);

    dJ=J0-J;

    if Direction=="Newton"
        p=p+gammaNewton*dp;
    elseif Direction=="gradient"
        p=p+gammaSD*g0SD;
    end

    pUpperViolation=pub-p;
    if any(pUpperViolation<0)
        fprintf("p violates pub! \n")
    end

    fprintf("%3i:(%s) \t gamma=%g \t J=%g \t J0=%g \t sub-obtimality=%g \t |dp|/|p|=%g \t J/J0=%g \n",iNewton,Direction,gamma,J,J0,SubOptimality,dpNorm,J/J0)

    Jvector(iNewton)=J0;
    Jvector(iNewton+1)=J;
    SlopeVector(iNewton)=slope0;
    gammaVector(iNewton)=gamma;
    SubOptimalityVector(iNewton)=SubOptimality;
    GradNormVector(iNewton)=norm(g0)/sqrt(numel(g0));


    drawnow

    % exit?
    if iNewton>(MaxNewtonSteps-1)
        fprintf("maximum number of iterations reached. \n")
        break
    end

    if SubOptimality<SubOptimalityTolerance
        fprintf("subtolerance reached. \n")
        break
    end

    if dJ<dJTolerance
        fprintf("J tolerance reached. \n")
        break
    end

    if dpNorm< dpTolerance
        fprintf("dp tolerance reached. \n")
        break
    end




end

I=~isnan(Jvector); Jvector=Jvector(I);
itVector=0:(numel(Jvector)-1) ; itVector=itVector(:);

I=~isnan(GradNormVector); GradNormVector=GradNormVector(I); GradNormVector=GradNormVector(:) ;  GradNormVector=[GradNormVector;NaN]; % Make sure it has the same length as itVector
I=~isnan(gammaVector); gammaVector=gammaVector(I); gammaVector=gammaVector(:) ;  gammaVector=[gammaVector;NaN]; % Make sure it has the same length as itVector

figIt=FindOrCreateFigure("J iteration") ; clf(figIt)
yyaxis left
semilogy(itVector,Jvector,"ob-",LineWidth=2) ;
ylabel("$J$",Interpreter="latex") ;
yyaxis right
semilogy(itVector,GradNormVector,"or-",LineWidth=2) ;
ylabel("$\|\nabla J \|$",Interpreter="latex") ;
xlabel("Iteration")

figItgamma=FindOrCreateFigure("gamma iteration") ; semilogy(itVector,gammaVector,"-k",MarkerFaceColor="r",MarkerEdgeColor="b",Marker="diamond")
ylabel("$\gamma$",Interpreter="latex") ;
xlabel("iteration")
title("step size in backtracking")

itRestart=max(RunInfo.Inverse.Iterations);

RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;itVector+itRestart];
RunInfo.Inverse.J=[RunInfo.Inverse.J;Jvector];
RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;GradNormVector];
RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gammaVector];




end

%%