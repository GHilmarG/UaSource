function [UserVar,RunInfo,LSF1,l,LSF1qx,LSF1qy,Residual]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)
%%
%
% $$\partial_t f +  \mathbf{v} \cdot \nabla f  - \nabla \cdot (\kappa \nabla f) = c \, \|(\nabla f)\|$$
%
%
%  df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f)
%
%    df/dt + (u-cx) df/dx + (v-cy) df/dy - div (kappa grad f) = 0
%
%
% Note: Tv, Lv and Pv are calculated only at the beginning of the NR iteration.
%

narginchk(7,8)
nargoutchk(4,7)

persistent iCalls


if ~CtrlVar.LevelSetMethod
    LSF1=F1.LSF;
    l=[];
    return
end



if isempty(iCalls)
    iCalls=0 ;
end
iCalls=iCalls+1;


if ~isempty(BCs.LSFL)
    L=BCs.LSFL;
    Lrhs=BCs.LSFrhs ;
else
    MLC=BCs2MLC(CtrlVar,MUA,BCs);
    L=MLC.LSFL ; Lrhs=MLC.LSFRhs ;
end


if nargin==7 || isempty(l) || (numel(l)~=numel(Lrhs))
    l=Lrhs*0 ;
end
dl=l*0 ;
dLSF=F1.LSF*0;
BCsError=0;

% make sure initial point is feasible
F1.LSF(BCs.LSFFixedNode)=BCs.LSFFixedValue;
F0.LSF(BCs.LSFFixedNode)=BCs.LSFFixedValue;

if CtrlVar.LSF.C
    if isempty(F1.LSFqx) || numel(F1.LSFqx)~=MUA.Nnodes
        F1.LSFqx=zeros(MUA.Nnodes,1);
        F1.LSFqy=zeros(MUA.Nnodes,1);
    end

    if isempty(F0.LSFqx)  || numel(F0.LSFqx)~=MUA.Nnodes
        F0.LSFqx=zeros(MUA.Nnodes,1);
        F0.LSFqy=zeros(MUA.Nnodes,1);
    end


    LSF1qx=F1.LSFqx;
    LSF1qy=F1.LSFqy;
else
    LSF1qx=[] ; LSF1qy=[] ;
end

iteration=0 ; rWork=inf ; rForce=inf; r=inf ; CtrlVar.NRitmin=0 ; gamma=1; rRatio=1;
RunInfo.LevelSet.SolverConverged=false; BackTrackInfo.iarm=0;

while true


    % exit criterion
    % First check if iteration>1, and backstep small and rRatio
    % small. There could be a good reason for this. For example we have
    % already reached the NR minimum and no further reduction is
    % possible due to rounding effect or becuase the cost function is
    % scaled in some ways that make further reduction irrelevant for
    % the solution.


    if rRatio<0.01 && rWork>1e-20 % OK I'm hardwiring in here an option of addional iteration
        % The argument being that we might most likely still be in the second-order convergence
        % and not limited by numerical rounding errors. (Consider taking this out once
        % all is fine.)

        ResidualsCriteria=false; % If there was a significant reduction in last iteration, do continue.

    elseif iteration>1 ...
            && (gamma < max(CtrlVar.LSFExitBackTrackingStepLength,CtrlVar.BacktrackingGammaMin) ...
            || rRatio>0.99)

        % This exist criterion applies if in last iteration either
        % backstep or the reduction was very small. This indicates
        % difficulties with convergence, but this might simply be
        % due to rounding errors making a further reduction
        % impossible.

        ResidualsCriteria=(rWork<CtrlVar.LSFDesiredWorkAndForceTolerances(1)  && rForce<CtrlVar.LSFDesiredWorkAndForceTolerances(2))...
            && (rWork<CtrlVar.LSFDesiredWorkOrForceTolerances(1)  || rForce<CtrlVar.LSFDesiredWorkOrForceTolerances(2))...
            && iteration >= CtrlVar.NRitmin;


    else

        % This is the otherwise default exit criterion
        ResidualsCriteria=(rWork<CtrlVar.LSFAcceptableWorkAndForceTolerances(1)  && rForce<CtrlVar.LSFAcceptableWorkAndForceTolerances(2))...
            && (rWork<CtrlVar.LSFAcceptableWorkOrForceTolerances(1)  || rForce<CtrlVar.LSFAcceptableWorkOrForceTolerances(2))...
            && iteration >= CtrlVar.NRitmin;


    end


    if iteration>=CtrlVar.LevelSetSolverMaxIterations && (rRatio>0.9)
        fprintf('LevelSetEquationNewtonRaphson: Maximum number of NR iterations (%i) reached. \n ',CtrlVar.LevelSetSolverMaxIterations)
        RunInfo.LevelSet.SolverConverged=false;
        break
    end


    if ResidualsCriteria
        fprintf('LevelSetEquationNewtonRaphson: NR iteration converged in %i iterations with rForce=%g and rWork=%g \n',iteration,rForce,rWork)
        RunInfo.LevelSet.SolverConverged=true;
        break
    end


    if RunInfo.BackTrack.Converged==0 || gamma<CtrlVar.LSFExitBackTrackingStepLength
        if CtrlVar.InfoLevelNonLinIt>=1
            fprintf('LevelSetEquationNewtonRaphson: Backtracking stagnated or step too small. Exiting after %i iterations with rWork=%g and rForce=%g \n',iteration,rWork,rForce)
        end
        RunInfo.LevelSet.SolverConverged=false;
        break
    end

    iteration=iteration+1 ;

    % [UserVar,R,K,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,F0.LSF,F0.c,F0.ub,F0.vb,F1.LSF,F1.c,F1.ub,F1.vb,F0.LSFqx,F0.LSFqy,LSF1qx,LSF1qy);
    [UserVar,R,K,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,F0,F1);


    if ~isempty(L)

        frhs=-R-L'*l        ;  % This needs to be identical to what is defined in the CalcCostFunctionLevelSetEquation
        grhs=Lrhs-L*F1.LSF  ;  % Here the argument is that frhs has the units: [\varphi] area/time
        % while grhs has the units [\varphi], where [\varphi] are the units of
        % the level-set function itself, ie [\varphi]
    else
        frhs=-R;
        grhs=[];
    end

    if issymmetric(K)
        [dLSF,dl]=solveKApeSymmetric(K,L,frhs,grhs,dLSF,dl,CtrlVar);
    else
        [dLSF,dl]=solveKApe(K,L,frhs,grhs,dLSF,dl,CtrlVar);
    end

    dLSF=full(dLSF);


    if CtrlVar.LSF.C  && ~isempty(Qx)  % consistent formulation

        if ~isfield(MUA,'dM') || isempty(MUA.dM)
            MUA.dM=decomposition(MUA.M,'chol','upper') ;
        end

        LSF1qx=MUA.dM\Qx ;
        LSF1qy=MUA.dM\Qy ;
        Residual=MUA.dM\Rv ;
    else
        LSF1qx=[] ; LSF1qy=[] ; Residual=[] ;
    end

    if any(isnan(dLSF))
        save TestSave
        error('LevelSetEquationNewtonRaphson:NaNinSolution','NaN in the solution for dLSF')
    end

    Func=@(gamma) CalcCostFunctionLevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dLSF,dl,BCs);

    gamma=0 ; [r0,~,~,rForce0,rWork0,D20]=Func(gamma);
    gamma=1 ; [r1,~,~,rForce1,rWork1,D21]=Func(gamma);

    slope0=-2*r0 ;

    if CtrlVar.LevelSetInfoLevel>=1000
        CtrlVar.InfoLevelBackTrack=1000 ;
    end


    CtrlVar.BacktrackingGammaMin=CtrlVar.LSFExitBackTrackingStepLength;
    [gamma,r,BackTrackInfo]=BackTracking(slope0,1,r0,r1,Func,CtrlVar);
    [r1Test,~,~,rForce,rWork,D2]=Func(gamma);



    if CtrlVar.LevelSetInfoLevel>=100 && CtrlVar.doplots==1
        nnn=30;
        gammaTestVector=zeros(nnn,1) ; rForceTestvector=zeros(nnn,1);  rWorkTestvector=zeros(nnn,1); rD2Testvector=zeros(nnn,1);
        Upper=2.2;
        Lower=-0.5 ;
        if gamma>0.7*Upper ; Upper=2*gamma; end
        parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:genericError');
        parfor I=1:nnn
            gammaTest=(Upper-Lower)*(I-1)/(nnn-1)+Lower
            [rTest,~,~,rForceTest,rWorkTest,D2Test]=Func(gammaTest);
            gammaTestVector(I)=gammaTest ; rForceTestvector(I)=rForceTest; rWorkTestvector(I)=rWorkTest;  rD2Testvector(I)=D2Test;
        end

        gammaZero=min(abs(gammaTestVector)) ;
        if gammaZero~=0
            [rTest,~,~,rForceTest,rWorkTest,D2Test]=Func(0);
            gammaTestVector(nnn+1)=0 ; rForceTestvector(nnn+1)=rForceTest; rWorkTestvector(nnn+1)=rWorkTest;  rD2Testvector(nnn+1)=D2Test;
        end

        [gammaTestVector,ind]=unique(gammaTestVector) ; rForceTestvector=rForceTestvector(ind) ; rWorkTestvector=rWorkTestvector(ind) ; rD2Testvector=rD2Testvector(ind) ;
        [gammaTestVector,ind]=sort(gammaTestVector) ; rForceTestvector=rForceTestvector(ind) ; rWorkTestvector=rWorkTestvector(ind) ; rD2Testvector=rD2Testvector(ind) ;


        SlopeForce=-2*rForce0;
        SlopeWork=-2*rWork0;
        SlopeD2=-D20;
        CtrlVar.MinimisationQuantity=CtrlVar.LSFMinimisationQuantity;
        [ForceFig,WorkFig]=PlotCostFunctionsVersusGamma(CtrlVar,RunInfo,gamma,r,iteration,"-LSF-",...
            gammaTestVector,rForceTestvector,rWorkTestvector,rD2Testvector,...
            SlopeForce,SlopeWork,SlopeD2,rForce,rWork,D2);

    end



    F1.LSF=F1.LSF+gamma*dLSF;
    l=l+gamma*dl;

    rRatio=r/r0;
    if CtrlVar.LevelSetInfoLevel>=1
        if ~isempty(L)
            BCsError=norm(Lrhs-L*F1.LSF);
        end

        fprintf(CtrlVar.fidlog,'Level-Set(TLP/%i%i%i):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , rForce=%-14.7g , rWork=%-14.7g , BCsError=%-14.7g \n ',...
            CtrlVar.LSF.T,CtrlVar.LSF.L,CtrlVar.LSF.P,...
            iteration,BackTrackInfo.iarm,gamma,rRatio,r0,r,rForce,rWork,BCsError);
    end

end

LSF1=F1.LSF ; % Because I don't return F1

RunInfo.LevelSet.iCount=RunInfo.LevelSet.iCount+1;

if numel(RunInfo.LevelSet.time) < RunInfo.LevelSet.iCount
    RunInfo.LevelSet.time=[RunInfo.LevelSet.time;RunInfo.LevelSet.time+NaN];
    RunInfo.LevelSet.Iterations=[RunInfo.LevelSet.Iterations;RunInfo.LevelSet.Iterations+NaN];
    RunInfo.LevelSet.Residual=[RunInfo.LevelSet.Residual;RunInfo.LevelSet.Residual+NaN];
    RunInfo.LevelSet.BackTrackSteps=[RunInfo.LevelSet.BackTrackSteps;RunInfo.LevelSet.BackTrackSteps+NaN];
    RunInfo.LevelSet.Phase=[RunInfo.LevelSet.Phase;strings(size(RunInfo.LevelSet.Phase))];
end



RunInfo.LevelSet.time(RunInfo.LevelSet.iCount)=CtrlVar.time;
RunInfo.LevelSet.Iterations(RunInfo.LevelSet.iCount)=iteration ;
RunInfo.LevelSet.Residual(RunInfo.LevelSet.iCount)=r;
RunInfo.LevelSet.BackTrackSteps( RunInfo.LevelSet.iCount)=BackTrackInfo.iarm ;
RunInfo.LevelSet.Phase(RunInfo.LevelSet.iCount)=CtrlVar.LevelSetPhase;



%     if CtrlVar.LevelSetInfoLevel>=1  && CtrlVar.LevelSetInfoLevel<10
%         if ~isempty(L)
%             BCsError=norm(Lrhs-L*F1.LSF);
%         end
%         fprintf(CtrlVar.fidlog,'Level-Set:%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , rForce=%-14.7g , rWork=%-14.7g , BCsError=%-14.7g \n ',...
%             iteration,BackTrackInfo.iarm,gamma,r/r0,r0,r,rForce,rWork,BCsError);
%     end
%



end


