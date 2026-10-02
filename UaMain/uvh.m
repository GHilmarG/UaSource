



function [UserVar,RunInfo,F1,l1,BCs1,dt]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)


narginchk(9,9)


if ~isfield(CtrlVar,"CminWidth") || ~isfield(CtrlVar,"AGlenminWidth") 
    CtrlVar=SetSmoothFloorWidths(CtrlVar,F0);
end

switch CtrlVar.uvhSolver

    case {"Root-Finding Newton","uvh_GaussNewton_fsolve"}


        [UserVar,RunInfo,F1,l1,BCs1,dt]=uvhRootFinding(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1);


    case "Least-Squares Gauss-Newton"

        [UserVar,RunInfo,F1,l1,BCs1,dt]=uvhLSQ(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1);
        RunInfo.Forward.uvhIterations(CtrlVar.CurrentRunStepNumber)=1;



    otherwise

        error("CaseFellThrough")

end


return

%%

CompareWithLSQ=false;

if CompareWithLSQ

    [UserVar,RunInfo,F1test,l1,BCs1,dt]=uvhLSQ(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1);

    UaPlots(CtrlVar,MUA,F1,"-uv-",FigureTitle="uv root (default)")
    UaPlots(CtrlVar,MUA,F1test,"-uv-",FigureTitle="uv lsq (MATLAB)")
    UaPlots(CtrlVar,MUA,F1test,[F1.ub-F1test.ub F1.vb-F1test.vb] ,FigureTitle="uv (diff)")


    UaPlots(CtrlVar,MUA,F1,F1.h,FigureTitle="h root (default)")
    UaPlots(CtrlVar,MUA,F1,F1.h,FigureTitle="h root (default)")
    UaPlots(CtrlVar,MUA,F1test,F1test.h,FigureTitle="h lsq (MATLAB)")
    UaPlots(CtrlVar,MUA,F1test,F1.h-F1test.h,FigureTitle="h diff")
end

end
