function [UserVar,RunInfo,h1,l]=MassContinuityEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)
 
 
%%
%
% $$\rho \partial h/\partial t + \partial  (\rho u h)/ \partial x + \partial  (\rho v h)/\partial y - \nabla \cdot (\kappa \nabla h ) = a(h)$$
%
%
%    
    narginchk(7,8)
    nargoutchk(6,6)

    
    MLC=BCs2MLC(CtrlVar,MUA,BCs);
    L=MLC.hL ; Lrhs=MLC.hRhs ;
    if nargin==7 || isempty(l) || (numel(l)~=numel(Lrhs))
        l=Lrhs*0 ; 
    end
    dl=l*0 ;
    dh=F1.h*0;
    BCsError=0;
    
    % make sure initial point is feasable
    F1.h(BCs.hFixedNode)=BCs.hFixedValue;
    F0.h(BCs.hFixedNode)=BCs.hFixedValue;
    
     
    iteration=0 ; rWork=inf ; rForce=inf; CtrlVar.NRitmin=0 ; gamma=1; rRatio=1;
    RunInfo.h.SolverConverged=false;
    
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
                && (gamma < max(CtrlVar.hExitBackTrackingStepLength,CtrlVar.BacktrackingGammaMin) ...
                || rRatio>0.99)
            
            % This exist criterion applies if in last iteration either
            % backstep or the reduction was very small. This indicates
            % difficulties with convergence, but this might simply be
            % due to rounding errors making a further reduction
            % impossible.
            
            ResidualsCriteria=(rWork<CtrlVar.hDesiredWorkAndForceTolerances(1)  && rForce<CtrlVar.hDesiredWorkAndForceTolerances(2))...
                && (rWork<CtrlVar.hDesiredWorkOrForceTolerances(1)  || rForce<CtrlVar.hDesiredWorkOrForceTolerances(2))...
                && iteration >= CtrlVar.NRitmin;
            
            
        else
            
            % This is the otherwise default exit criterion
            ResidualsCriteria=(rWork<CtrlVar.hAcceptableWorkAndForceTolerances(1)  && rForce<CtrlVar.hAcceptableWorkAndForceTolerances(2))...
                && (rWork<CtrlVar.hAcceptableWorkOrForceTolerances(1)  || rForce<CtrlVar.hAcceptableWorkOrForceTolerances(2))...
                && iteration >= CtrlVar.NRitmin;
            
            
        end
        
        
        if iteration>=CtrlVar.hSolverMaxIterations && (rRatio>0.9) 
            fprintf('hEquationNewtonRaphson: Maximum number of NR iterations (%i) reached. \n ',CtrlVar.hSolverMaxIterations)
            RunInfo.h.SolverConverged=false;
            break
        end
        
        
        if ResidualsCriteria 
            fprintf('hEquationNewtonRaphson: NR iteration converged in %i iterations with rForce=%g and rWork=%g \n',iteration,rForce,rWork)
            RunInfo.h.SolverConverged=true;
            break
        end
        
        
        if RunInfo.BackTrack.Converged==0 || gamma<CtrlVar.hExitBackTrackingStepLength 
            if CtrlVar.InfoLevelNonLinIt>=1
                fprintf('hEquationNewtonRaphson: Backtracking stagnated or step too small. Exiting after %i iterations with rWork=%g and rForce=%g \n',iteration,rWork,rForce)
            end
            RunInfo.h.SolverConverged=false;
            break
        end
        
        iteration=iteration+1 ;
        
        [UserVar,R,K]=MassContinuityEquationAssembly(UserVar,CtrlVar,MUA,F0.h,F0.rho,F0.ub,F0.vb,F0.a,F0.dadh,F1.h,F1.c,F1.ub,F1.vb,F1.a,F1.dadh);
        
        if ~isempty(L)
   
            frhs=-R-L'*l        ;  % This needs to be identical to what is defined in the CalcCostFunctionhEquation
            grhs=Lrhs-L*F1.h; % Here the argument is that frhs has the units: [\varphi] area/time
                                % while grhs has the units [\varphi], where [\varphi] are the units of 
                                % the level-set function itself, ie [\varphi]
        else
            frhs=-R;
            grhs=[];
        end
        
        [dh,dl]=solveKApe(K,L,frhs,grhs,dh,dl,CtrlVar);
        dh=full(dh);

        
        if any(isnan(dh))
            save TestSave
            error('hEquationNewtonRaphson:NaNinSolution','NaN in the solution for dh')
        end
        
        Func=@(gamma) CalcCostFunctionhEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dh,dl);

        gamma=0 ; [r0,~,~,rForce0,rWork0,D20]=Func(gamma);
        gamma=1 ; [r1,~,~,rForce1,rWork1,D21]=Func(gamma);
        
        slope0=-2*r0 ;
        
        if CtrlVar.hInfoLevel>=1000
            CtrlVar.InfoLevelBackTrack=1000 ;
        end
        
       
        CtrlVar.BacktrackingGammaMin=CtrlVar.hExitBackTrackingStepLength; 
        [gamma,r,BackTrackInfo]=BackTracking(slope0,1,r0,r1,Func,CtrlVar);
        [r1Test,~,~,rForce,rWork,D2]=Func(gamma);
        
        
        % Testing
        if gamma<1e-5 
           fprintf(" hm \n ") 
        end
        
        if CtrlVar.hInfoLevel>=100 && CtrlVar.doplots==1
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
            CtrlVar.MinimisationQuantity=CtrlVar.hMinimisationQuantity;
            [ForceFig,WorkFig]=PlotCostFunctionsVersusGamma(CtrlVar,RunInfo,gamma,r,iteration,"-h-",...
                gammaTestVector,rForceTestvector,rWorkTestvector,rD2Testvector,...
                SlopeForce,SlopeWork,SlopeD2,rForce,rWork,D2);
            
        end
        
        
        
        F1.h=F1.h+gamma*dh;
        l=l+gamma*dl;
        
        rRatio=r/r0; 
        if CtrlVar.hInfoLevel>=1
            if ~isempty(L)
                BCsError=norm(Lrhs-L*F1.h);
            end
           
            fprintf(CtrlVar.fidlog,'Level-Set(TLP/%i%i%i):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , rForce=%-14.7g , rWork=%-14.7g , BCsError=%-14.7g \n ',...
                 CtrlVar.h.T,CtrlVar.h.L,CtrlVar.h.P,...
                iteration,BackTrackInfo.iarm,gamma,rRatio,r0,r,rForce,rWork,BCsError);
        end
        
    end
    
    h1=F1.h ; % Because I don't return F1
    
    RunInfo.h.iCount=RunInfo.h.iCount+1;
    
    if numel(RunInfo.h.time) < RunInfo.h.iCount
       RunInfo.h.time=[RunInfo.h.time;RunInfo.h.time+NaN];
       RunInfo.h.Iterations=[RunInfo.h.Iterations;RunInfo.h.Iterations+NaN];
       RunInfo.h.Residual=[RunInfo.h.Residual;RunInfo.h.Residual+NaN];
       RunInfo.h.BackTrackSteps=[RunInfo.h.BackTrackSteps;RunInfo.h.BackTrackSteps+NaN];
       RunInfo.h.Phase=[RunInfo.h.Phase;strings(size(RunInfo.h.Phase))]; 
    end
    
    
    
    RunInfo.h.time(RunInfo.h.iCount)=CtrlVar.time;   
    RunInfo.h.Iterations(RunInfo.h.iCount)=iteration ; 
    RunInfo.h.Residual(RunInfo.h.iCount)=r;
    RunInfo.h.BackTrackSteps( RunInfo.h.iCount)=BackTrackInfo.iarm ; 
    RunInfo.h.Phase(RunInfo.h.iCount)=CtrlVar.hPhase;
    
    
    
    
end

