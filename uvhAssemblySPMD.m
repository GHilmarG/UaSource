function [UserVar,RunInfo,R,K]=uvhAssemblySPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1)




if isempty(CtrlVar.Parallel.uvhAssembly.spmd.nWorkers)
    poolobj = gcp('nocreate');
    CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=poolobj.NumWorkers;
end

nW=CtrlVar.Parallel.uvhAssembly.spmd.nWorkers;

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
    
    % all ele based quantities need to be transferred accross the the lab
    M{labindex}.connectivity=MUA.connectivity(Partition{labindex},:);
    M{labindex}.Nele=numel(Partition{labindex});
    M{labindex}.Deriv=MUA.Deriv(Partition{labindex},:,:,:);
    M{labindex}.DetJ=MUA.DetJ(Partition{labindex},:);
    M{labindex}.EleAreas=MUA.EleAreas(Partition{labindex},:);
     
    [~,~,rr,kk]=uvhMatrixAssembly(UserVar,RunInfo,CtrlVar,M{labindex},F0,F1);
end

% I don't see how this can be done differently, but this is slow
R=rr{1};
for iWorker=2:nW
    R=R+rr{iWorker};
end

if ~isempty(kk)
    K=kk{1};
    for iWorker=2:nW
        K=K+kk{iWorker};
    end
end

end