

%%
%
% 2023-01: General conclusions.
%
% For full matrices it is easy to get significant speed up my simply creating distributed arrays. 
%
% For sparse matrices it is also easy to get significant speedup, provided that matrix has density of 0.1 or more
%
% For matrices that arise with Ua the density is however much smaller, or about 0.0001, and no significant speedup is
% observed. 
%
%
%%

% delete(gcp('nocreate')); parpool('Processes',8)

% delete(gcp('nocreate')); parpool('Threads',8)

% distributed not supported in thread-based parallel pool (F2024a)

% https://uk.mathworks.com/help/parallel-computing/benchmarking-a-b.html

%% full solver
%
% Very simple approach to evaluate parallel performance when solving full systems.
%
% This clearly does speed things up, although exact performance gain will depend on various factors.
%

nWVector=[1 2 4 6 8 10 12];
tSolveVector=nWVector+nan;
perfVector=nWVector+nan;
iCount=0;
N=20000;
A = rand( N );
b = sum(A,2) ;


for nW=nWVector

    delete(gcp('nocreate'));
    parpool('Processes',nW) ;

    ADist=distributed(A);
    bDist=distributed(b); 

    tSolve=tic ; % Start of timed region
    x = ADist \ bDist; % Solve the linear system
    tSolve = toc(tSolve) ; % End of timed region
    perf = ( 2/3 * N^3 + 3/2 * N^2 ) / tSolve / 1.e9 ;

    iCount=iCount+1;
    tSolveVector(iCount)=tSolve;
    perfVector(iCount)=perf; 
    fprintf("nW=%i \t time %f sec \t GFlops %f \n ",nW,tSolve,perf)


end

figure(30)
yyaxis left
plot(nWVector,tSolveVector,"o-b")
ylabel("time for solve (sec)",Interpreter="latex")
yyaxis right
plot(nWVector,perfVector,"o-r")
ylabel("Gigaflops",Interpreter="latex")

xlabel("\# workers",Interpreter="latex") ; 
title("direct solver, full system, 20k$\times$20k ",Interpreter="latex")





%% sparse solver

%% full solver
%
% Very simple approach to evaluate parallel performance when solving full systems.
%
% This clearly does speed things up, although exact performance gain will depend on various factors.
%

filename=[]; 
filename="solveKApePIGTWGuvh250896time0k19NRit2";


nWVector=[1 2 4 6 8 10 12];
nWVector=[1  4  8  10 ] ;
tSolveVector=nWVector+nan;
perfVector=nWVector+nan;
iCount=0;


if isempty(filename)

    n=20;
    N=1000*n;
    density=0.1 ;
    A = sprand( N , N , density);
    b = sum(A,2) ;

else

    load(filename,"A")
    b = sum(A,2) ;
    N=size(A,1) ;
    density=nnz(A)/N^2; 

end

for nW=nWVector

    delete(gcp('nocreate'));
    parpool('Processes',nW) ;

    ADist=distributed(A);
    bDist=distributed(b); 

    tSolve=tic ; % Start of timed region
    x = ADist \ bDist; % Solve the linear system
    tSolve = toc(tSolve) ; % End of timed region
    perf = ( 2/3 * N^3 + 3/2 * N^2 ) / tSolve / 1.e9 ;

    iCount=iCount+1;
    tSolveVector(iCount)=tSolve;
    perfVector(iCount)=perf; 
    fprintf("nW=%i \t time %f sec \t GFlops %f \n ",nW,tSolve,perf)


end

figure(500)
% yyaxis left
plot(nWVector,tSolveVector,"o-b")
ylabel("time for solve (sec)",Interpreter="latex")
% yyaxis right ; plot(nWVector,perfVector,"o-r") ; ylabel("Gigaflops",Interpreter="latex")
xlabel("\# workers",Interpreter="latex") ; 


title(sprintf("direct solver, sparse system, %ik $\\times$ %ik, density=%g ",round(N/10e3),round(N/10e3),density),Interpreter="latex")






%%
