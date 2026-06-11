function [Ruv,Kuv,Tint,Fext]=uvAssemblyParfeval(CtrlVar,MUA,F)


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
%spmd (0,nW)

tic
for iWorker=1:nW
    % Build M directly on the workers to avoid communication
    
        
    M{iWorker}=MUA;
    M{iWorker}.connectivity=MUA.connectivity(Partition{iWorker},:);
    M{iWorker}.Nele=numel(Partition{iWorker});
    M{iWorker}.Deriv=MUA.Deriv(Partition{iWorker},:,:,:);
    M{iWorker}.DetJ=MUA.DetJ(Partition{iWorker},:);

end
toc




   
tic    
for iWorker=1:nW

    Future(iWorker)=parfeval(@uvMatrixAssembly,2,CtrlVar,M{iWorker},F);

end
toc


tic
[Ruv,Kuv] = fetchOutputs(Future(1));

for iWorker=2:nW
    [rrTemp,kkTemp] = fetchOutputs(Future(iWorker));
    Ruv=Ruv+rrTemp ;
    Kuv=Kuv+kkTemp ;
end
toc

Tint=[] ; Fext=[]; 





end