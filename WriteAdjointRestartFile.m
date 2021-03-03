function WriteAdjointRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo,InvStartValues,Priors,Meas,BCsAdjoint,InvFinalValues)


fprintf(CtrlVar.fidlog,'Saving inverse-run restart file: %s \n ',CtrlVar.Inverse.NameOfRestartOutputFile);

CtrlVarInRestartFile=CtrlVar;
UserVarInRestartFile=UserVar;
time=CtrlVar.time;
dt=CtrlVar.dt;
save(CtrlVar.Inverse.NameOfRestartOutputFile,...
    'CtrlVarInRestartFile','UserVarInRestartFile','MUA','BCs','F','GF','l','RunInfo',...
    'InvStartValues','Priors','Meas','BCsAdjoint','InvFinalValues','time','dt','-v7.3');


F.x=MUA.coordinates(:,1) ; 
F.y=MUA.coordinates(:,2) ; 

if CtrlVar.AGlenisElementBased
    xA=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1));
    yA=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2));
else
    xA=MUA.coordinates(:,1);
    yA=MUA.coordinates(:,2);
end

if CtrlVar.CisElementBased
    xC=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1));
    yC=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2));
else
    xC=MUA.coordinates(:,1);
    yC=MUA.coordinates(:,2);
end

   
if CtrlVar.Inverse.SaveSlipperinessEstimateInSeperateFile
    fprintf(CtrlVar.fidlog,' saving sliding-law parameters  in file %s \n ',CtrlVar.NameOfFileForSavingSlipperinessEstimate) ;
    C=F.C;
    m=F.m;
    muk=F.muk ; 
    q=F.q ; 
    save(CtrlVar.NameOfFileForSavingSlipperinessEstimate,'C','m','q','muk','xC','yC','MUA','CtrlVarInRestartFile')
end

if CtrlVar.Inverse.SaveAGlenEstimateInSeperateFile
    fprintf(CtrlVar.fidlog,' saving AGlen and m in file %s \n ',CtrlVar.NameOfFileForSavingAGlenEstimate) ;
    AGlen=F.AGlen;
    n=F.n;
    save(CtrlVar.NameOfFileForSavingAGlenEstimate,'AGlen','n','xA','yA','MUA','CtrlVarInRestartFile')
end

end