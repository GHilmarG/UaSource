
%%


NumWorkers=8 ;

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

load("uvAtilde.mat","Atilde","btilde")

spparms ('spumoni', 0)
Gain=nan(32,1);


for NumWorkers=[1 2 3 4  6 8 10 12]

    ParPool = gcp('nocreate') ;

    if isempty(ParPool)

        parpool('Processes',NumWorkers);

    elseif (ParPool.NumWorkers~=NumWorkers)

        delete(gcp('nocreate'));
        parpool('Processes',NumWorkers);

    end

    for iRepeat=1:3

        tSeq=tic;

        x=Atilde\btilde;

        tSeq=toc(tSeq) ;

        tDist=tic;

        AtildeDist=distributed(Atilde);
        btildeDist=distributed(btilde);
        x=AtildeDist\btildeDist;
        x=gather(x);
        tDist=toc(tDist);



        fprintf("%i : tSeq=%g \t tDist=%g \t gain=%g \n",NumWorkers,tSeq,tDist,tSeq/tDist)
    end

    Gain(NumWorkers)=tSeq/tDist;
end

% rather odd behaviour: about factor 2 gain irrespectivly of number of workers...
% even just one worker results in about gain of 2...Something else must be going on here.
figure(100) ; plot(Gain,'or')


%%



%%
load("solveKApePIGTWGuvh250896.mat","A","B","CtrlVar","f","g","x0","y0")


CtrlVar.Parallel.Distribute=false ;

tSeq=tic;
[x,y]=solveKApe(A,B,f,g,x0,y0,CtrlVar);
tSeq=toc(tSeq);

CtrlVar.Parallel.Distribute=true ;

tDist=tic;
[x,y]=solveKApe(A,B,f,g,x0,y0,CtrlVar);
tDist=toc(tDist);


% hm, seems better to do the distribution just ahead of the \ operation, ie within ABfgPreEliminate
tDist2=tic;
Adist=distributed(A); Bdist=distributed(B);
[x,y]=solveKApe(Adist,Bdist,f,g,x0,y0,CtrlVar);
tDist2=toc(tDist2);

fprintf("%i : tSeq=%g \t tDist=%g \t tDist2=%g \t gain=%g \n",NumWorkers,tSeq,tDist,tDist2, tSeq/tDist)

%%




