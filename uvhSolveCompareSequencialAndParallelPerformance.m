   
function [UserVar,RunInfo,F,l,BCs,dt]=uvhSolveCompareSequencialAndParallelPerformance(UserVar,RunInfo,CtrlVar,MUA,F0,F,l,BCs)

%% First do the uvh solve using whichever parallel performance options the user has already set.

CtrlVar.Parallel.isTest=false ; % set this to false to suppress further performance comparisons within the uvh call.
RunInfo.CPU.Assembly.uvh=0 ;  RunInfo.CPU.Solution.uvh=0 ; % Reset this cumulative sum of CPU sec used for assembly and linsolve.

CtrlVar.Parallel.uvhAssembly.spmd.isOn=true; CtrlVar.Parallel.uvAssembly.spmd.isOn=true; CtrlVar.Parallel.Distribute=true;
MUA=UpdateMUA(CtrlVar,MUA) ;

tParallel=tic ;
[UserVar,RunInfo]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F,l,l,BCs);  
tParallel=toc(tParallel);

tAssemblyPar=RunInfo.CPU.Assembly.uvh ; 
tSolvePar=RunInfo.CPU.Solution.uvh ; 

%% And now for comparison turn all parallel options off, and repeat the uvh solve.
CtrlVar.Parallel.uvhAssembly.spmd.isOn=false; CtrlVar.Parallel.uvAssembly.spmd.isOn=false; CtrlVar.Parallel.Distribute=false;

RunInfo.CPU.Assembly.uvh=0 ;  RunInfo.CPU.Solution.uvh=0 ; % Again, teset this cumulative sum of CPU sec used for assembly and linsolve.


tSeq=tic ;

[UserVar,RunInfo,F,l,BCs,dt]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F,l,l,BCs);
tSeq=toc(tSeq);

tAssemblySeq=RunInfo.CPU.Assembly.uvh ; 
tSolveSeq=RunInfo.CPU.Solution.uvh ; 

%% Summarize info

fprintf("\n \n Comparision between sequential and parallel uvh solve using the user-defined parallel options set in DefineInitialInputs.m \n ")
fprintf(" uvh Solve:  \t tSeq=%10f sec \t tPar=%10f sec \t tSeq/tPar=%5f \n",tSeq,tParallel,tSeq/tParallel)
fprintf("  assembly:  \t tSeq=%10f sec \t tPar=%10f sec \t tSeq/tPar=%5f \n",tAssemblySeq,tAssemblyPar,tAssemblySeq/tAssemblyPar)
fprintf("  linsolve:  \t tSeq=%10f sec \t tPar=%10f sec \t tSeq/tPar=%5f \n \n",tSolveSeq,tSolvePar,tSolveSeq/tSolvePar)

end 