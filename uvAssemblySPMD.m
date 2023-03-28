function [Ruv,Kuv,Tint,Fext]=uvAssemblySPMD(CtrlVar,MUA,F)


if isempty(CtrlVar.Parallel.uvAssembly.spmd.nWorkers)
    poolobj = gcp;
    CtrlVar.Parallel.uvAssembly.spmd.nWorkers=poolobj.NumWorkers;
end

nW=CtrlVar.Parallel.uvAssembly.spmd.nWorkers;

Partition=cell(nW,1);
M=cell(nW,1);

%% create element lists for each partition
% Use round to make sure that there are exactly nW partitions
% and ensure that all elements are included.
N=round(MUA.Nele/nW) ; i1=1 ; i2=N;
for iWorker=1:(nW-1)
    Partition{iWorker}=i1:i2 ;
    i1=i2+1 ;
    i2=i2+N ;
end

i2=MUA.Nele;
Partition{nW}=i1:i2 ;
%

MUA.dM=[] ; 
spmd (0,nW) 
    % Build M directly on the workers to avoid communication
    
        
    M{spmdIndex}=MUA;
    M{spmdIndex}.connectivity=MUA.connectivity(Partition{spmdIndex},:);
    M{spmdIndex}.Nele=numel(Partition{spmdIndex});
    M{spmdIndex}.Deriv=MUA.Deriv(Partition{spmdIndex},:,:,:);
    M{spmdIndex}.DetJ=MUA.DetJ(Partition{spmdIndex},:);
    [rr,kk]=uvMatrixAssembly(CtrlVar,M{spmdIndex},F);
    rrsum = spmdPlus(rr,1);
    kksum = spmdPlus(kk,1);
end


Ruv=rrsum{1};
Kuv=kksum{1};

Tint=[] ; Fext=[]; 



end