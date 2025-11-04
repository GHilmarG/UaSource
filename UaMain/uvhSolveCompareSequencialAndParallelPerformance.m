   
function [UserVar,RunInfo,F,l,BCs,dt]=uvhSolveCompareSequencialAndParallelPerformance(UserVar,RunInfo,CtrlVar,MUA,F0,F,l,BCs)



ParPool = gcp('nocreate');  % check if parpool exists, but do not create one if it does not exist already


if ~isempty(ParPool)

    CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=ParPool.NumWorkers;
    
    parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:genericError');
    parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:SaveNotSupported');

    if  CtrlVar.Parallel.uvhAssembly.spmd.isOn
        CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=ParPool.NumWorkers;
    end

    if CtrlVar.Parallel.uvAssembly.spmd.isOn  
        CtrlVar.Parallel.uvAssembly.spmd.nWorkers=ParPool.NumWorkers;
    end

end




%% First do the uvh solve using whichever parallel performance options the user has already set.


% make a copy of user settings
User.CtrlVar.Parallel.uvhAssembly.spmd.isOn=CtrlVar.Parallel.uvhAssembly.spmd.isOn;
User.CtrlVar.Parallel.uvAssembly.spmd.isOn=CtrlVar.Parallel.uvAssembly.spmd.isOn;
User.CtrlVar.Parallel.Distribute=CtrlVar.Parallel.Distribute;

CtrlVar.Parallel.isTest=false ; % set this to false to suppress further performance comparisons within the uvh call.
RunInfo.CPU.Assembly.uvh=0 ;  RunInfo.CPU.Solution.uvh=0 ; % Reset this cumulative sum of CPU sec used for assembly and linsolve.

%CtrlVar.Parallel.uvhAssembly.spmd.isOn=true; CtrlVar.Parallel.uvAssembly.spmd.isOn=true; CtrlVar.Parallel.Distribute=true;
tic
MUA=UpdateMUA(CtrlVar,MUA) ;
toc

tParallel=tic ;
[UserVar,RunInfo]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F,l,l,BCs);  
tParallel=toc(tParallel);

tAssemblyPar=RunInfo.CPU.Assembly.uvh ; 
tSolvePar=RunInfo.CPU.Solution.uvh ; 

%% And now for comparison turn all parallel options off, and repeat the uvh solve.
CtrlVar.Parallel.uvhAssembly.spmd.isOn=false; CtrlVar.Parallel.uvAssembly.spmd.isOn=false; CtrlVar.Parallel.Distribute=false;

RunInfo.CPU.Assembly.uvh=0 ;  RunInfo.CPU.Solution.uvh=0 ; % Again, reset this cumulative sum of CPU sec used for assembly and linsolve.


tSeq=tic ;
[UserVar,RunInfo,F,l,BCs,dt]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F,l,l,BCs);
tSeq=toc(tSeq);

tAssemblySeq=RunInfo.CPU.Assembly.uvh ; 
tSolveSeq=RunInfo.CPU.Solution.uvh ; 

%% Summarize info
[status,hostname]=system('hostname');
hostname=strtrim(convertCharsToStrings(hostname)) ;
fprintf("\n \n ------------- Comparison between sequential and parallel uvh solve using the user-defined parallel options set in DefineInitialInputs.m ------------ \n ")
fprintf("Machine: %s :  #Elements=%i \t #Nodes=%i \t #Workers=%i \t \n",hostname,MUA.Nele,MUA.Nnodes,CtrlVar.Parallel.uvhAssembly.spmd.nWorkers)
fprintf(" Note: Parallel options not switched on by user are not used in the parallel solve.\n")
fprintf(" uvh-Solve (total time):  \t tSeq=%10f sec \t tPar=%10f sec \t tSeq/tPar=%5f \n",tSeq,tParallel,tSeq/tParallel)

if User.CtrlVar.Parallel.uvhAssembly.spmd.isOn
    fprintf("              assembly:  \t tSeq=%10f sec \t tPar=%10f sec \t tSeq/tPar=%5f \n",tAssemblySeq,tAssemblyPar,tAssemblySeq/tAssemblyPar)
end

if User.CtrlVar.Parallel.Distribute
    fprintf("  linsolve:  \t tSeq=%10f sec \t tPar=%10f sec \t tSeq/tPar=%5f \n \n",tSolveSeq,tSolvePar,tSolveSeq/tSolvePar)
else
   fprintf("  distributed option was not used when solving the matrix equation \n ")
end

end 