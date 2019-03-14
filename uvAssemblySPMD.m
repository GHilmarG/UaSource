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

spmd (0,nW) 
    % Build M directly on the workers to avoid communication
    
        
    M{labindex}=MUA;
    M{labindex}.connectivity=MUA.connectivity(Partition{labindex},:);
    M{labindex}.Nele=numel(Partition{labindex});
    M{labindex}.Deriv=MUA.Deriv(Partition{labindex},:,:,:);
    M{labindex}.DetJ=MUA.DetJ(Partition{labindex},:);
    [rr,kk]=uvMatrixAssembly(CtrlVar,M{labindex},F);
    rrsum = gplus(rr,1);
    kksum = gplus(kk,1);
end


Ruv=rrsum{1};
Kuv=kksum{1};

Tint=[] ; Fext=[]; 



end