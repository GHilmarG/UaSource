



function [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F,BCs)

narginchk(4,4)

persistent iCount

if isempty(iCount)
    iCount=0;
end


if isempty(CtrlVar.Parallel.uvAssembly.spmd.nWorkers)
    poolobj = gcp;
    CtrlVar.Parallel.uvAssembly.spmd.nWorkers=poolobj.NumWorkers;
end

nW=CtrlVar.Parallel.uvAssembly.spmd.nWorkers;


% outside of spmd  M is  composite
% inside of spmd M is struct on each worker

MUA.dM=[] ;


tAssembly=tic;
MUAworkers=MUA.workers; 
spmd (nW)
    [rr,kk]=uvMatrixAssemblySSTREAM(CtrlVar,MUAworkers,F,BCs);
end
tAssembly=toc(tAssembly);



tSum=tic ;
spmd (nW)
    rrsum = spmdPlus(rr,1);
    kksum = spmdPlus(kk,1);
end
% rrsum and kksum are composites
% Ruv and Kuv are double sparse
Ruv=rrsum{1}; Kuv=kksum{1};

tSum=toc(tSum) ;


% Delete on workers. This should not be needed, but for some reason the permanence in threaded environment degrades with the
% number of calls. This is an attempt to reset this, but this did not have the desired effect.  This performance degrade does
% not happen in process environment.  (11 Jan 2024)
% spmd ; rr=[] ; kk=[] ; rrsum=[] ; kksum=[] ; end


Tint=[] ; Fext=[];


if CtrlVar.Parallel.isTest  && CtrlVar.InfoLevelCPU>=10

    iCount=iCount+1;
    fprintf("uvMatrixAssemblySSTREAM_SPMD (%i): SPMD Assembly %f sec. \t Summing up results from workers %f sec.  \n",...
        iCount,tAssembly,tSum)
end


end