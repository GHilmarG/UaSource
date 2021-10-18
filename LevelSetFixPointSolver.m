
function [UserVar,RunInfo,LSF,l]=LevelSetFixPointSolver(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)



%clearvars ; load TestSaveFixPoint.mat



fprintf("\n\n\n\n\n\n")

FixPointConverged=0 ; 
CtrlVar.LevelSetFixPointSolverApproach="-initial fix point (IFP) then pseudo time stepping (PTS) followed by final fix point (FFP)-";

RunInfo.LevelSet.SolverConverged=false ;
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


if ~FixPointConverged && contains(CtrlVar.LevelSetFixPointSolverApproach,"PTS") 

    % pseudo time stepping
    l=[] ;
    CtrlVar.LSF.P=1 ; CtrlVar.LSF.T=1 ; CtrlVar.LSF.L=0 ;
    CtrlVar.LevelSetTheta=1;

    dtOld=CtrlVar.dt ;

    N=0 ; Factor=5 ; TotalTime=0 ; dLSFdt=inf ;
    Nmax=50;
    dLSFdtMax= 1 ;  TotalTimeMax=100000;

    while N<=Nmax  && dLSFdt > dLSFdtMax  && TotalTime<TotalTimeMax

        N=N+1;


        fprintf("\n \n =======================   Pseudo-time steping #%i with dt=%f and Total Time=%f with Res=%f ====================================== \n ",N,CtrlVar.dt,TotalTime,dLSFdt)

        [UserVar,RunInfo,LSF,l]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);


        if RunInfo.LevelSet.SolverConverged
            TotalTime=TotalTime+CtrlVar.dt ;
            dLSFdt=max(F1.LSF-LSF)/CtrlVar.dt ;
            F1.LSF=LSF ; F0.LSF=F1.LSF ;
            if RunInfo.LevelSet.Iterations(RunInfo.LevelSet.iCount)<5
                CtrlVar.dt=Factor*CtrlVar.dt;
            end
        else
            CtrlVar.dt=CtrlVar.dt/Factor;
            dLSFdt=inf ;
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

