%


%


TestCase="-direct-" ;
% TestCase="-compare-" ;
% TestCase="-best-"  ;


switch TestCase


    case "-direct-"

        % Comaprision with Direct solver

        load solveKApePIGTWGuvh250896.mat
        tic
        [x,y,tolA,tolB]=ABfgPreEliminate(CtrlVar,A,B,f,g) ;
        toc
        x0=x;
        y0=y;


    case "-compare-"

        % Iterative method comparision


        Peq=[] ;
        Req=[] ;
        Ceq=[] ;
        x0=[];
        y0=[];

        CtrlVar.InfoLevelLinSolve=100;
        [x,y,tolA,tolB,Peq,Req,Ceq]=ABfgPreEliminateIterativeMethodComparision(CtrlVar,A,B,f,g,x0,y0,Peq,Req,Ceq) ;


    case "-best-"

%         From the method comparsion it is concluded that ilutp+dissect+gmrs is the best approach.
%         Using equlibriate does not improve the convergence signficantly, and takes lot of time (30 sec, for the typical example
%         used here).
% 
%         So fastes approach appearse to be:
% 
%         1) Preconditoner based in ilutp with droptolerance around 1e-6 or so
%         2) dissect  (a must to limit memory)
%         3) gmres with tol=1e-15 and maxit around 30 (Usually only five or so needed)


        L=[]; U=[] ; x0=[] ; y0=[] ; perm=[] ; 
        



        % load solveKApePIGTWGuvh250896.mat ;
        load("solveKApePIGTWGuvh250896time0k19NRit2.mat","A","B","CtrlVar","f","g","x0","y0")

        % A=distributed(A) ; B=distributed(B) ;  f=distributed(f) ; g=distributed(g) ; x0=distributed(x0) ; y0=distributed(y0) ;  % this does not work because dissect does not support distributed arrays
        % A=gpuArray(A) ; B=gpuArray(B) ;  f=gpuArray(f) ; g=gpuArray(g) ; x0=gpuArray(x0) ; y0=gpuArray(y0) ;  % this does not work because dissect does not support distributed arrays

        CtrlVar.InfoLevelLinSolve=100;

        x0=[] ; y0=[] ; 
        fprintf("\n\n\n-----------------------------------------------------------------------------------------\n\n")
        tstart1=tic ;
        [x,y,tolA,tolB,L,U,perm]=ABfgPreEliminateIterative(CtrlVar,A,B,f,g,x0,y0,L,U,perm) ;
        tend1=toc(tstart1) ; 

        
        
        
        
        CtrlVar.InfoLevelLinSolve=100;

        fprintf("\n\n\n-----------------------------------------------------------------------------------------\n\n")
        % x0=x+1e-8*abs(x).*rand(length(x0),1) ; % for this slight modification of x0, a repeated solve is fast and does not require
        % new LU factorisation.
        
        load("solveKApePIGTWGuvh250896time0k19NRit3.mat","A","B","CtrlVar","f","g","x0","y0")
        % A=gpuArray(A) ; B=gpuArray(B) ;  f=gpuArray(f) ; g=gpuArray(g) ; x0=gpuArray(x0) ; y0=gpuArray(y0) ;
        CtrlVar.InfoLevelLinSolve=100;
        % x0=x ; y0=y; % L=[] ; U=[] ; perm=[] ; 
        tstart2=tic ;
        [x,y,tolA,tolB,L,U,perm]=ABfgPreEliminateIterative(CtrlVar,A,B,f,g,x0,y0,L,U,perm) ;
        tend2=toc(tstart2) ; 


        fprintf("tend1=%f \t tend2=%f \n",tend1,tend2)


    otherwise

        fprintf(" Which case? \n")

end



