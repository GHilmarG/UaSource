




function MUAworkers=BuildMuaWorkers(CtrlVar,MUA,MUAworkers)


if ~CtrlVar.Parallel.BuildWorkers  % Building the worker data can be suppressed: AdaptMesh.m does so for the intermediate meshes.
                                   % (31 Oct 2026) The default is now true; UpdateMUA.m only calls this function when needed.

    
    return

end



poolobj = gcp('nocreate') ;   % (31 Oct 2026) never create a pool here (gcp would): without a pool, no worker data are built
if isempty(poolobj)
    return
end
if isempty(CtrlVar.Parallel.uvhAssembly.spmd.nWorkers)
    CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=poolobj.NumWorkers;
end

nW=CtrlVar.Parallel.uvhAssembly.spmd.nWorkers;

%%





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

% (26 Oct 2026) The worker data are always built here: the callers decide when this is necessary. UpdateMUA.m keeps the key of
% the mesh and the number of workers for which the worker data were built (MUA.workersKey, see MuaWorkersKey.m), and checks with
% all(exist(MUA.workers)) that the Composite belongs to the current pool. The previous check here (number of nodes, and the element
% partitions compared on the workers) required communication with the workers, and did not detect a change of the mesh that leaves
% the numbers of nodes and elements unchanged. The input MUAworkers is only returned unchanged if CtrlVar.Parallel.BuildWorkers=false.

BuildWorkers=true ;

if BuildWorkers
    MUA.workers=[];   % Make sure to get rid of this composite, within the MUA structure ahead of the spmd call
    % First of all this field is no longer needed as it has be judged to be in an incorrect state and needs to be recreated.
    % Second, keeping it means that this composite, ie MUA.workers, is used within the body of the spmd statement and this is not
    % allowed in matlab and will create a warning.
    spmd (0,nW)

        % Build M directly on the workers to avoid communication

        MUAworkers.nod=MUA.nod;
        MUAworkers.nip=MUA.nip;

        MUAworkers.Nnodes=MUA.Nnodes;
        MUAworkers.points=MUA.points;
        MUAworkers.weights=MUA.weights;

        MUAworkers.coordinates=MUA.coordinates;
        %
        MUAworkers.Partition=Partition{spmdIndex} ;
        MUAworkers.connectivity=MUA.connectivity(Partition{spmdIndex},:);
        MUAworkers.Nele=numel(Partition{spmdIndex});
        MUAworkers.Deriv=MUA.Deriv(Partition{spmdIndex},:,:,:);
        MUAworkers.DetJ=MUA.DetJ(Partition{spmdIndex},:);
        MUAworkers.EleAreas=MUA.EleAreas(Partition{spmdIndex}) ;
        [MUAworkers.uvAssemblyPattern,MUAworkers.uvhAssemblyPattern]=AssemblyPatternCache(CtrlVar,MUAworkers) ;   % (17 Oct 2026) the uv pattern as well, for uvMatrixAssemblySSTREAM_SPMD.m   
                                                                                                                  % (16 Oct 2026) sparsity pattern of the sub-domain, for the faster pattern-based assembly


    end
end



end