function [UserVar,RunInfo,R,K,Tint,Fext]=uvhAssemblySPMD(UserVar,RunInfo,CtrlVar,MUA,F0,F1,ZeroFields)


nW=CtrlVar.NumWorkers;
N=ceil(MUA.Nele/nW) ;


i1=1 ;  i2=N;
for iWorker=1:nW
    
    Partition=[i1:i2] ;
    i1=i2+1 ;
    i2=min([i2+N,MUA.Nele]);
    
    M{iWorker}=MUA;
    M{iWorker}.connectivity=MUA.connectivity(Partition,:);
    M{iWorker}.Nele=numel(Partition);
    M{iWorker}.Deriv=MUA.Deriv(Partition,:,:,:);
    M{iWorker}.DetJ=MUA.DetJ(Partition,:);
    
end

Tint=[] ; Fext=[];

spmd
    [~,~,rr,kk]=uvhAssembly(UserVar,RunInfo,CtrlVar,M{labindex},F0,F1,ZeroFields);
end

R=rr{1}; K=kk{1};
for iWorker=2:nW
    R=R+rr{iWorker};
    K=K+kk{iWorker};
end


end