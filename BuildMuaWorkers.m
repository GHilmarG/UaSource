



function MUAworkers=BuildMuaWorkers(CtrlVar,MUA,MUAworkers)


if isempty(CtrlVar.Parallel.uvhAssembly.spmd.nWorkers)
    poolobj = gcp ;
    CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=poolobj.NumWorkers;
end

nW=CtrlVar.Parallel.uvhAssembly.spmd.nWorkers;


%% create element lists for each partition
% Use round to make sure that there are exactly nW partitions
% and ensure that all elements are included.

Partition=cell(nW,1);
N=round(MUA.Nele/nW) ; i1=1 ; i2=N;
for iWorker=1:(nW-1)
    Partition{iWorker}=i1:i2 ;
    i1=i2+1 ;
    i2=i2+N ;
end

i2=MUA.Nele;
Partition{nW}=i1:i2 ;

% outside of spmd  M is  composite
% inside of spmd M is struct on each worker

MUA.dM=[] ;


if isempty(MUAworkers) || numel(MUAworkers)==0
    spmd (0,nW)

        % Build M directly on the workers to avoid communication

        MUAworkers.nod=MUA.nod;
        MUAworkers.nip=MUA.nip;

        MUAworkers.Nnodes=MUA.Nnodes;
        MUAworkers.points=MUA.points;
        MUAworkers.weights=MUA.weights;

        MUAworkers.coordinates=MUA.coordinates;
        %
        MUAworkers.connectivity=MUA.connectivity(Partition{spmdIndex},:);
        MUAworkers.Nele=numel(Partition{spmdIndex});
        MUAworkers.Deriv=MUA.Deriv(Partition{spmdIndex},:,:,:);
        MUAworkers.DetJ=MUA.DetJ(Partition{spmdIndex},:);
        MUAworkers.EleAreas=MUA.EleAreas(Partition{spmdIndex}) ;

    end
end



end