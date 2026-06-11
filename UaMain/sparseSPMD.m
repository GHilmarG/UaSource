function K=sparseSPMD(CtrlVar,iK,jK,Kval,ni,nj,nW)

%%
%
% This is a simple attempt to use SPMD to create sparse arrays on workers.
%
% This seems to be able to speed things in a threaded pool, at least if the array is large
%
%
%%

n=numel(iK);

tSPMD=tic;

Partition=cell(nW,1);
N=round(n/nW) ;
i1=1 ; i2=N;

for iWorker=1:(nW-1)
    Partition{iWorker}=i1:i2 ;
    i1=i2+1 ;
    i2=i2+N ;
end

i2=n;
Partition{nW}=i1:i2 ;


tSparse=tic;
% Building the co-distributed arrays directly on the workers. I expect that there are fancier ways of doing this. This is
% basically a manually created co-distributor, but I can't understand MATLAB documentation on how to do this otherwise. 
spmd (nW)

    Kworkers=sparse(iK(Partition{spmdIndex}),jK(Partition{spmdIndex}),Kval(Partition{spmdIndex}),ni,nj);

end
tSparse=toc(tSparse);


tSum=tic;
spmd (0,nW)
    KworkersSum = spmdPlus(Kworkers,1);
end
K=KworkersSum{1};
tSum=toc(tSum) ;
tSPMD=toc(tSPMD);

if ~isempty(CtrlVar)
    if CtrlVar.Parallel.isTest

        tSeq=tic;
        Kseq=sparse(iK,jK,Kval,ni,nj); 
        tSec=toc(tSeq) ;
     
        fprintf("sparseSPMD:  normest(Kseq)-KSPMD)/normest(diag(Kseq))=%g \n",normest(Kseq-K)/normest(Kseq))
        fprintf("sparseSPMD:  SPMD sparse on workers %f sec. \t Summing up results from workers %f sec.  \n",tSparse,tSum)
        fprintf("sparseSPMD:  sparse outside SPMD in %f sec, \t sparse SPMD in %f sec \t Speedup=%f \n ",tSec,tSPMD,tSec/tSPMD)


    end
end

end



