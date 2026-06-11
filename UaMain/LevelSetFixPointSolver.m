
function [UserVar,RunInfo,LSF,l]=LevelSetFixPointSolver(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)



%clearvars ; load TestSaveFixPoint.mat



fprintf("\n LSF: fix point solver \n")

FixPointConverged=0 ;



RunInfo.LevelSet.SolverConverged=false ;
lIFP=[] ;
if contains(CtrlVar.LevelSetFixPointSolverApproach,"IFP")
    % Fix point
    l=[] ;
    CtrlVar.LSF.P=1 ; CtrlVar.LSF.T=0 ; CtrlVar.LSF.L=0 ;
    CtrlVar.LevelSetTheta=1; CtrlVar.LevelSetSolverMaxIterations=20;
    [UserVar,RunInfo,LSF,l]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
    if RunInfo.LevelSet.SolverConverged
        F1.LSF=LSF ;
        lIFP=l ;   FixPointConverged=1 ;
    else
        lIFP=[] ;   FixPointConverged=0 ;
    end

end


if ~FixPointConverged || contains(CtrlVar.LevelSetFixPointSolverApproach,"PTS")

    % pseudo time stepping
    l=[] ;
    CtrlVar.LSF.P=1 ; CtrlVar.LSF.T=1 ; CtrlVar.LSF.L=0 ;
    CtrlVar.LevelSetTheta=1;

    dtOld=CtrlVar.dt ;

    N=0 ;
    FactorUp=1.5 ; FactorDown=2 ; iUpLast=true ;
    TotalTime=0 ; dLSFdt=inf ;



    Nmax=CtrlVar.LevelSetPseudoFixPointSolverMaxIterations ;
    dLSFdtMax= CtrlVar.LevelSetPseudoFixPointSolverTolerance ;
    TotalTimeMax=1e8*CtrlVar.dt;

    while true

        N=N+1;


        fprintf("\n \n =======================   Pseudo-time steping #%i with dt=%f and Total Time=%f with Res=%f ====================================== \n ",N,CtrlVar.dt,TotalTime,dLSFdt)

        [UserVar,RunInfo,LSF,l]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);


        if RunInfo.LevelSet.SolverConverged
            TotalTime=TotalTime+CtrlVar.dt ;
            dLSFdt=max(F1.LSF-LSF)/CtrlVar.dt ;
            F1.LSF=LSF ; F0.LSF=F1.LSF ;

            if  dLSFdt < dLSFdtMax

                fprintf('Exiting pseudo-step fix-point LSF initialisation because desired tolerance on dLSF/dt has been met.\n')
                fprintf('Desired tolerance was set to %f, and on exit tolerance is %f. \n',dLSFdtMax,dLSFdt)
                break

            end


            if RunInfo.LevelSet.Iterations(RunInfo.LevelSet.iCount)<4

                if iUpLast
                    iUpLast=false;
                else
                    CtrlVar.dt=FactorUp*CtrlVar.dt;
                    iUpLast=true;
                end
            end

        else
            CtrlVar.dt=CtrlVar.dt/FactorDown;
            dLSFdt=inf ;
        end


        if  N >= Nmax

            fprintf('Exiting pseudo-step fix-point LSF initialisation because maximum number of iterations (%i) reached.\n',Nmax)
            break

        end

        if  TotalTime>TotalTimeMax

            fprintf('Exiting pseudo-step fix-point LSF initialisation because total integration time (%f) reached.\n',TotalTimeMax)
            break

        end



    end


    fprintf("\n \n =========== Pseudo-time steping exit after #%i iterations, with dt=%f and Total Time=%f with Res=%f ====================================== \n ",N,CtrlVar.dt,TotalTime,dLSFdt)
    CtrlVar.dt=dtOld;
end




if ~FixPointConverged && contains(CtrlVar.LevelSetFixPointSolverApproach,"FFP")
    fprintf("\n\n\n Now do the fix point \n\n\n")
    l=lIFP;
    CtrlVar.LSF.P=1 ; CtrlVar.LSF.T=0 ; CtrlVar.LSF.L=0 ;
    CtrlVar.LevelSetTheta=1; CtrlVar.LevelSetSolverMaxIterations=20;
    [UserVar,RunInfo,LSF,l]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
end


end

