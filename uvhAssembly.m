function [UserVar,RunInfo,R,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1)
    
    narginchk(6,6)

    nargoutchk(4,4)
    tAssembly=tic;
    
    if CtrlVar.Parallel.uvhAssembly.spmd.isOn
        [UserVar,RunInfo,R,K]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
    else
        [UserVar,RunInfo,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
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



    %%

    if CtrlVar.Parallel.isTest && ~CtrlVar.uvhMatrixAssembly.Ronly

        tSeq=tic;
        [UserVar,RunInfo,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
        tSeq=toc(tSeq) ;

        CtrlVar.Parallel.uvhAssembly.spmd.isOn=true ;
        MUA=UpdateMUA(CtrlVar,MUA);  % just in case workers have not been set up
        tSPMD=tic;
        [UserVar,RunInfo,Rspmd,Kspmd]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
        tSPMD=toc(tSPMD);


        if isempty(CtrlVar.Parallel.uvhAssembly.spmd.nWorkers)
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=0;
            else
                CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=poolobj.NumWorkers;
            end
        end

        fprintf('\n ----------------------------- Info on parallel uvh SPMD assembly performance : Nele=%i  \t nWorkers=%i \n',MUA.Nele,CtrlVar.Parallel.uvAssembly.spmd.nWorkers)
        fprintf('#Ele=%i \t SPMD used for uvh assembly:  tSeq=%f \t tSPMD=%f \t speedup=%g \n',MUA.Nele,tSeq,tSPMD,tSeq/tSPMD) ;
        fprintf(' R-Rspmd=%g \t K-Kspmd=%g   \n',full(norm(R-Rspmd)/norm(R)),normest(K-Kspmd)/normest(K))
        fprintf(' ----------------------------- \n')




    end
    %

end
