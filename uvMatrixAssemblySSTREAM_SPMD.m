



function [Ruv,Kuv,Tint,Fext,MUAworkers]=uvMatrixAssemblySSTREAM_SPMD(CtrlVar,MUA,F,MUAworkers)

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






%% create element lists for each partition
% Use round to make sure that there are exactly nW partitions
% and ensure that all elements are included.
tPartition=tic;
Partition=cell(nW,1);
N=round(MUA.Nele/nW) ; i1=1 ; i2=N;
for iWorker=1:(nW-1)
    Partition{iWorker}=i1:i2 ;
    i1=i2+1 ;
    i2=i2+N ;
end

i2=MUA.Nele;
Partition{nW}=i1:i2 ;
tPartition=toc(tPartition);
%


% outside of spmd  M is  composite
% inside of spmd M is struct on each worker

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
% rrsum and kksum are composites
% Ruv and Kuv are double sparse
Ruv=rrsum{1}; Kuv=kksum{1};

tSum=toc(tSum) ;


% Delete on workers. This should not be needed, but for some reason the permanence in threaded environment degrades with the
% number of calls. This is an attempt to reset this, but this did not have the desired effect.  This performance degrade does
% not happen in process environment.  (11 Jan 2024)
spmd ; rr=[] ; kk=[] ; rrsum=[] ; kksum=[] ; end


Tint=[] ; Fext=[];


if CtrlVar.Parallel.isTest

    iCount=iCount+1;
    fprintf("uvMatrixAssemblySSTREAM_SPMD (%i): Creating partition arrays %f sec. \t Building arrays on workers %f sec. \t SPMD Assembly %f sec. \t Summing up results from workers %f sec.  \n",...
        iCount,tPartition,tBuild,tAssembly,tSum)
end


end