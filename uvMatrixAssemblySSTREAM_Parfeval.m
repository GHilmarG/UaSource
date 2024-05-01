



function [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM_Parfeval(CtrlVar,MUA,F,BCs)

%%
%
% [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM_Parfeval(CtrlVar,MUA,F)
%
% SSTREAM matrix assembly done by using parfeval over element partitions (domains)
%
%
%%


if isempty(CtrlVar.Parallel.uvAssembly.spmd.nWorkers)
    poolobj = gcp;
    CtrlVar.Parallel.uvAssembly.spmd.nWorkers=poolobj.NumWorkers;
end

nW=CtrlVar.Parallel.uvAssembly.spmd.nWorkers;


%% create element lists for each partition
Partition=cell(nW,1);
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
%%
MUA.dM=[]; 

f(1:nW)=parallel.FevalFuture;

tParfeval=tic;
for idx=1:nW

    f(idx)=parfeval(backgroundPool,@uvMatrixAssemblySSTREAMpartitionTriplets,6,CtrlVar,MUA,F,Partition{idx}); 
end
tParfeval=toc(tParfeval);


tFetching=tic ;
[iK,jK,Kval,iR,Tval,Fval]=fetchOutputs(f); %  Each output from f is concatenated, so here I get 7 vectors with all triplet values from all the workers
tFetching=toc(tFetching);

cancel(f)  ; % Not sure if this is needed

tSparse=tic ;
dof=2;   neq=dof*MUA.Nnodes;


Kuv=sparseUA(iK,jK,Kval,neq,neq);
% Kuv=sparseSPMD(CtrlVar,iK,jK,Kval,neq,neq,nW) ;

Tint=sparseUA(iR,1,Tval);
Fext=sparseUA(iR,1,Fval);
Ruv=Tint-Fext;
tSparse=toc(tSparse);

if CtrlVar.Parallel.isTest
    fprintf("uvAssemblyParfevalSSTREAM: parfeval call %f sec. Fetching outputs from workers %f sec. \t Creating sparse matrices %f sec.  \n",tParfeval,tFetching,tSparse)
end





end