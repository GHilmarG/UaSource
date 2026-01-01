

function [p,UserVar,RunInfo]=BruteForceHessianInversion(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub)


%%
%load("TestSaveH.mat","func","p0","CtrlVar","iRange","MUA","F")  ;


%%

p=p0;

nNewtonSteps=1;
Jvector=nan(nNewtonSteps+1,1);
SlopeVector=nan(nNewtonSteps+1,1);
gammaVector=nan(nNewtonSteps+1,1);
SubOptimalityVector=nan(nNewtonSteps+1,1);
GradNormVector=nan(nNewtonSteps+1,1);

nPar=numel(p);

iRange=1:nPar; % all of them

CtrlVar.LineSearchAllowedToUseExtrapolation=true;


SubOptimalityTolerance=1e-7;
dJTolerance=0.0;
dpTolerance=0.0;

iNewton=0;
lmin=0.1; lEnd=0; 
while true

    iNewton=iNewton+1;

    [Hsparse,Hfull,g0,J0] = CalcBruteForceHessian(func,p,CtrlVar,iRange) ;

    if isnan(J0)
        error("BruteForceHessianInversion:J0IsNaN","NaN in J0")
    end


    lStart=max(lEnd,lmin);
    [Hfull,lEnd]=CheckIfHessianIsSPDandIfNotMakeItSo(Hfull,MUA,lStart) ;

    dp=Hfull\(-g0);
 

    if anynan(dp)
        fprintf("Solving the Newton system resulted in nan. \n")
        error("BruteForceHessianInversion:dpIsNaN","NaN in dp")
    end

    if anynan(p)
        fprintf("p contains nan. \n")
        error("BruteForceHessianInversion:pIsNaN","NaN in p")
    end


    % Since the Hessian may have been modified, it is not clear what a sensible first step could be
    gamma=1;
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
          CtrlVar.NewtonAcceptRatio=0.1 ;CtrlVar.BacktrackingGammaMin=0.001;
          [gammaNewton,JNewton,BackTrackInfo]=BackTracking(slope0,gamma,J0,J1,Func,CtrlVar);
    else
          fprintf("Slope at origin in Newton direction is positive! \n")
          gammaNewton=nan;
          JNewton=inf;
    end



    % I also do Steepest decent, and then compare.
    % This is because calculating the Hessian is so expensive and, in comparison, the line-search cheap

    % Then maybe modify H and solve until slope negative, or just go for the gradient direction
    nNodes=size(MUA.M,1);
    nH=size(Hfull,1);
    if nH==nNodes
        EYE=MUA.M ;
    elseif nH==2*nNodes
        EYE=blkdiag(MUA.M,MUA.M) ;
    else
        error("wrong dimentions")
    end

    g0SD=EYE*g0; % pre-multiplying
    Func=@(gamma) func(p-gamma*g0SD); % switching to the direction of the (negative) gradient
    slope0=-g0'*g0SD;
    gamma=-0.1*J0/slope0;
    J1=Func(gamma);
  
    CtrlVar.NewtonAcceptRatio=0.1 ;CtrlVar.BacktrackingGammaMin=gamma/100;
    [gammaSD,JSD,BackTrackInfo]=BackTracking(slope0,gamma,J0,J1,Func,CtrlVar);

    fprintf("====> JNewton/J0=%g \t JSD/J=%g \n",JNewton/J0,JSD/J0)

    if JNewton  < JSD
        fprintf("Newton step resulted in  greater reduction than steepest decent \n ")
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



    PlotCostVersusStepSizeAlongNewtonDirection(func,p,dp,g0,gammaNewton,JNewton,g0SD,gammaSD,JSD);

    dJ=J0-J;

    if Direction=="Newton"
        p=p+gammaNewton*dp;
    elseif Direction=="gradient"
        p=p-gammaSD*g0SD;
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
    if iNewton>(nNewtonSteps-1)
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


itRestart=max(RunInfo.Inverse.Iterations); 

RunInfo.Inverse.Iterations=[RunInfo.Inverse.Iterations;itVector+itRestart];
RunInfo.Inverse.J=[RunInfo.Inverse.J;Jvector];
RunInfo.Inverse.GradNorm=[RunInfo.Inverse.GradNorm;GradNormVector];
RunInfo.Inverse.StepSize=[RunInfo.Inverse.StepSize;gammaVector];





end

%%