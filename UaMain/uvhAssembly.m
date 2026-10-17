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

    if CtrlVar.InfoLevelCPU>=10
        if CtrlVar.Parallel.uvhAssembly.spmd.isOn
            if CtrlVar.uvhMatrixAssembly.Ronly
                fprintf("SPMD Force vector assembly in %g sec \n",tAssembly)
            else
                fprintf("SPMD Matrix assembly in %g sec \n",tAssembly)
            end
        else
            if CtrlVar.uvhMatrixAssembly.Ronly
                fprintf("Force vector assembly in %g sec \n",tAssembly)
            else
                fprintf("Matrix assembly in %g sec \n",tAssembly)
            end
        end
    end



 
end
