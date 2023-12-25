%


%%
%
% These are various ideas that I've tested while trying to get an interative solver performance that is better then direct
% solver
%
% I have never been able to get the iterative solver to be even close to the direct one in terms of performance.
%
% The best option so far was use of gmres with ilu preconditioner combined with dissect. This was 'only' about 2 to 5 times slover than
% the direct solver. Also found the performance of the iterative solver to be highly problem dependent.
%
% Seems impossible to speed this up using paralle options as disect not supported for distributed arrays (2023b)
% 
%
%%


TestCase="-direct-" ;
TestCase="-compare-" ;
TestCase="-best-"  ;

load("solveKApePIGTWGuvh250896.mat","A","B","CtrlVar","f","g","x0","y0")

switch TestCase


    case "-direct-"

        % Comaprision with Direct solver

        
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
        L1=[];
        U1=[];
        L1eq=[];
        U1eq=[];

        

        load("solveKApePIGTWGuvh250896time0k19NRit2.mat","A","B","CtrlVar","f","g")
        CtrlVar.InfoLevelLinSolve=100;
        
        
        % A=distributed(A) ; B=distributed(B) ; f=distributed(f) ; g=distributed(g); % equilibrate and dissect do not work with
        % either distributed nor gpuarrays

        x0=[];
        y0=[];

                        
        fprintf("\n\n\n-----------------------------------------------------------------------------------------\n\n")

        [x,y,tolA,tolB,Peq,Req,Ceq,L1,U1,L1eq,U1eq]=ABfgPreEliminateIterativeMethodComparision(CtrlVar,A,B,f,g,x0,y0,Peq,Req,Ceq,L1,U1,L1eq,U1eq) ; 

        % Now call again using previous preconditioners and equilibrated matrix

         fprintf("\n\n\n-----------------------------------------------------------------------------------------\n\n")

         % Now load a matrix from next NR iteration
         load("solveKApePIGTWGuvh250896time0k19NRit3.mat","A","B","CtrlVar","f","g") ; 
         CtrlVar.InfoLevelLinSolve=100;
         
         x0=x ; y0=y ; % use previous solution as an initial guess, and the previous preconditioners as well


        [x,y,tolA,tolB,Peq,Req,Ceq,L1,U1,L1eq,U1eq]=ABfgPreEliminateIterativeMethodComparision(CtrlVar,A,B,f,g,x0,y0,Peq,Req,Ceq,L1,U1,L1eq,U1eq) ; 


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
        x0=x ; y0=y; % L=[] ; U=[] ; perm=[] ; 
        tstart2=tic ;
        [x,y,tolA,tolB,L,U,perm]=ABfgPreEliminateIterative(CtrlVar,A,B,f,g,x0,y0,L,U,perm) ;
        tend2=toc(tstart2) ; 


        fprintf("tend1=%f sec \t tend2=%f sec \n",tend1,tend2)


    otherwise

        fprintf(" Which case? \n")

end



