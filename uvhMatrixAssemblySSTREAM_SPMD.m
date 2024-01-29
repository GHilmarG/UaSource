

function [UserVar,RunInfo,Ruvh,Kuvh]=uvhMatrixAssemblySSTREAM_SPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1)

narginchk(6,6)

if isempty(CtrlVar.Parallel.uvhAssembly.spmd.nWorkers)
    poolobj = gcp;
    CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=poolobj.NumWorkers;
end


nW=CtrlVar.Parallel.uvhAssembly.spmd.nWorkers;



MUAworkers=MUA.workers;

spmd (nW) 

    [~,~,rr,kk]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,MUAworkers,F0,F1);

end

spmd (0,nW)
    rrsum = spmdPlus(rr,1);
    kksum = spmdPlus(kk,1);
end

% rrsum and kksum are composites
% Ruvh and Kuvh are double sparse

Ruvh=rrsum{1}; 
Kuvh=kksum{1};



end