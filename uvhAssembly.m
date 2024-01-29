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

    if CtrlVar.Parallel.isTest && ~CtrlVar.uvhMatrixAssembly.Ronly

        tSeq=tic;
        [UserVar,RunInfo,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
        tSeq=toc(tSeq) ;

        tSPMD1=tic;
        [UserVar,RunInfo,R1,K1]=uvhAssemblySPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1);
        tSPMD1=toc(tSPMD1);

        tSPMD2=tic;
        [UserVar,RunInfo,R2,K2]=uvhAssemblySPMD2(UserVar,RunInfo,CtrlVar,MUA,F0,F1);
        tSPMD2=toc(tSPMD2);


        fprintf('\n \n ----------------------------- \n')
        fprintf(' tSeq=%f \t tSPMD1=%f \t tSPMD2=%f \n',tSeq,tSPMD1,tSPMD2) ;

    
        fprintf(' R-R1=%g \t K-K1=%g   \n',full(norm(R-R1)/norm(R)),norm(diag(K)-diag(K1))/norm(diag(K))); 
        fprintf(' R-R2=%g \t K-K2=%g   \n \n',full(norm(R-R2)/norm(R)),norm(diag(K)-diag(K2))/norm(diag(K))); 
    


    end
    %

end
