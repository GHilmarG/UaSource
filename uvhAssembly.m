function [UserVar,RunInfo,R,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1)
    
    narginchk(6,6)

    nargoutchk(4,4)
    tAssembly=tic;
    
    if CtrlVar.Parallel.uvhAssembly.spmd.isOn
        [UserVar,RunInfo,R,K]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
    else
        [UserVar,RunInfo,R,K]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;
    end


    RunInfo.CPU.Assembly.uvh=RunInfo.CPU.Assembly.uvh+toc(tAssembly); 

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


        fprintf('\n ----------------------------- Info on parallel SPMD assembly performance : \n')
        fprintf('SPMD used for uvh assembly:  tSeq=%f \t tSPMD=%f \t tSeq/rSPMD=%f \n',tSeq,tSPMD,tSeq/tSPMD) ;
        fprintf(' R-Rspmd=%g \t K-Kspmd=%g   \n',full(norm(R-Rspmd)/norm(R)),norm(diag(K)-diag(Kspmd))/norm(diag(K))); 
        fprintf(' ----------------------------- \n')
        
    


    end
    %

end
