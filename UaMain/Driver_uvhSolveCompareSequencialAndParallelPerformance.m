

%%
load("../UaTests/Calving/PIG-TWG/uvhSolve-PIG-TWG-10km.mat","UserVar","CtrlVar","RunInfo","F","MUA","F0","l","BCs")

%%
delete(gcp('nocreate'))
parpool("processes",8)

MUA.workers=[] ; 
tic
MUA=UpdateMUA(CtrlVar,MUA) ;
toc

CtrlVar.Parallel.uvhAssembly.spmd.isOn=true;
CtrlVar.Parallel.uvAssembly.spmd.isOn=true; 
CtrlVar.Parallel.Distribute=true; 


[UserVar,RunInfo,F,l,BCs,dt]=uvhSolveCompareSequencialAndParallelPerformance(UserVar,RunInfo,CtrlVar,MUA,F0,F,l,BCs) ;

%%


