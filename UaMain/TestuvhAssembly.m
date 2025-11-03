delete(gcp('nocreate'))
poolobj=parpool('local',[1 2]);

%%
poolobj = gcp('nocreate');


%clc
fprintf('\n\n\n\n')
%load uvhAssembly.mat
load TestSaveuvhAssembly

%%
fprintf('\n\n--------------------------------------------------------------------------------------------------\n\n')
CtrlVar.Parallel.uvhAssembly.parfor=0;


tic
[UserVar,RunInfo,R,K,Tint,Fext]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,ZeroFields);
toc
fprintf(' sum(R)=%g \t trace(K)=%g \n',full(sum(R)),full(trace(K)))

% The spmd option has potential

CtrlVar.NumWorkers=poolobj.NumWorkers;
fprintf(' Number of workers=%i \n ',CtrlVar.NumWorkers)
tic
[UserVar,RunInfo,R,K,Tint,Fext]=uvhAssemblySPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1,ZeroFields);
toc
fprintf(' sum(R)=%g \t trace(K)=%g \n',full(sum(R)),full(trace(K)))

%% The parfor option is not good and apparantly always slower


nW=poolobj.NumWorkers;
N=ceil(MUA.Nele/nW) ; 

i1=1 ;  i2=N;

for I=1:nW

    Partition{I}=[i1:i2] ;
    i1=i2+1 ;
    i2=min([i2+N,MUA.Nele]);
end

% using parfor
neq=3*MUA.Nnodes;
R=sparseUA(neq,1);
K=sparseUA(neq,neq);

tic
parfor iWorker=1:numel(Partition)
    
    M{iWorker}=MUA;
    M{iWorker}.connectivity=MUA.connectivity(Partition{iWorker},:);
    M{iWorker}.Nele=numel(Partition{iWorker});
    M{iWorker}.Deriv=MUA.Deriv(Partition{iWorker},:,:,:);
    M{iWorker}.DetJ=MUA.DetJ(Partition{iWorker},:);
    [~,~,r,k,Tint,Fext]=uvhAssembly(UserVar,RunInfo,CtrlVar,M{iWorker},F0,F1,ZeroFields);
    R=R+r ; 
    K=K+k;
    
end
toc





