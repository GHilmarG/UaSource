function [Ruv,Kuv,Tint,Fext,MUAworkers]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F,MUAworkers)

narginchk(4,4)


if isempty(CtrlVar.Parallel.uvAssembly.spmd.nWorkers)
    poolobj = gcp;
    CtrlVar.Parallel.uvAssembly.spmd.nWorkers=poolobj.NumWorkers;
end

nW=CtrlVar.Parallel.uvAssembly.spmd.nWorkers;

Partition=cell(nW,1);



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


% outside of spmd  M is  composite
% inside of spmd M is struct

MUA.dM=[] ;



tBuild=tic;
if isempty(MUAworkers)
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

    end
end
tBuild=toc(tBuild) ;



tAssembly=tic;
spmd (0,nW)
    [rr,kk]=uvMatrixAssemblySSTREAM(CtrlVar,MUAworkers,F);
end
tAssembly=toc(tAssembly);



tSum=tic ;
spmd (0,nW)
    rrsum = spmdPlus(rr,1);
    kksum = spmdPlus(kk,1);
end
tSum=toc(tSum) ;


Ruv=rrsum{1};
Kuv=kksum{1};

Tint=[] ; Fext=[];


if CtrlVar.Parallel.isTest
    fprintf("uvMatrixAssemblySSTREAM_SPMD: Building arrays on workers %f sec. \t SPMD Assembly %f sec. \t Summing up results from workers %f sec.  \n",tBuild,tAssembly,tSum)
end


end