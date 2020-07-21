function [UserVar,RunInfo,R,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1)
    
    narginchk(6,6)
    nargoutchk(4,4)
    tAssembly=tic;
    
    if CtrlVar.Parallel.uvhAssembly.spmd.isOn
        [UserVar,RunInfo,R,K]=uvhAssemblySPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1);
    else
        [UserVar,RunInfo,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
    end
    
    
    RunInfo.CPU.Assembly.uvh=toc(tAssembly)+RunInfo.CPU.Assembly.uvh;
    
    %%
    
    if CtrlVar.Parallel.isTest
        
        tSeq=tic;
        [UserVar,RunInfo,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
        tSeq=toc(tSeq) ;
        
        tSPMD=tic;
        [UserVar,RunInfo,R2,K2]=uvhAssemblySPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1);
        tSPMD=toc(tSPMD);
        
        if ~CtrlVar.uvhMatrixAssembly.Ronly
            
            fprintf('\n \n ----------------------------- \n')
            fprintf(' tSeq=%f \t tSPMD=%f \t MUA.Nele=%i \n',tSeq,tSPMD,MUA.Nele) ;
            if CtrlVar.WriteRunInfoFile
                fprintf(RunInfo.File.fid,' tSeq=%f \t tSPMD=%f \t MUA.Nele=%i \t %f  \t %f \n',tSeq,tSPMD,MUA.Nele,full(norm(R-R2)/norm(R)),full(sum(abs(max(K)-max(K2)))/sum(abs(max(K))))) ;
            end
        end
        
    end
    %
    
end
