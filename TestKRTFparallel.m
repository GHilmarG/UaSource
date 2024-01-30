
%%


NumWorkers=12 ;

ParPool = gcp('nocreate') ;

if isempty(ParPool)

    parpool('Processes',NumWorkers)

elseif (ParPool.NumWorkers~=NumWorkers)

    delete(gcp('nocreate'))
    parpool('Processes',NumWorkers)

end

%%
N=10000;
A=rand(N,N);
B=rand(N,N);

%%
clc
tic
AA=distributed(A);
BB=distributed(B);
toc

tic
C=TestAddM(A,B);
t=toc;
trace(C)

tic

    CC=TestAddM(AA,BB);

tspmd=toc;
trace(CC)

fprintf('t=%f \t tspmd=%f \n',t,tspmd)







%%



load("TestSolveKApe.mat","A","B","f","g","x0","y0","CtrlVar")


%%

tSeq=tic;
[x,y]=ABfgPreEliminate(CtrlVar,A,B,f,g);
tSeq=toc(tSeq);

tDistribute=tic;
Adist=distributed(A);
Bdist=distributed(B);
fdist=distributed(f);
gdist=distributed(g);
tDistribute=toc(tDistribute);

tPar=tic;
[xDist,yDist]=ABfgPreEliminate(CtrlVar,Adist,Bdist,fdist,gdist);
x2=gather(xDist) ; y2=gather(yDist) ;   
tPar=toc(tPar);

fprintf("tSeq=%g \t tDistribute=%g \t tPar=%g \t gain=%g \n",tSeq,tDistribute,tPar,tSeq/(tDistribute+tPar))

% tSeq=26.2609 	 tDistribute=1.36496 	 tPar=13.7741 	 gain=1.73465  C23000099     8 workers
%%



