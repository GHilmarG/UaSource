function [UserVar,RunInfo,R,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)
    
    narginchk(8,8)

    nargoutchk(3,4)
    tAssembly=tic;
    
    if CtrlVar.Parallel.uvhAssembly.spmd.isOn
        [UserVar,RunInfo,R,K]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;
    else
        [UserVar,RunInfo,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;
    end

    tAssembly=toc(tAssembly);
    RunInfo.CPU.Assembly.uvh=RunInfo.CPU.Assembly.uvh+tAssembly;



 
end
