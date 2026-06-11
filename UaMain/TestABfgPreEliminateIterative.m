%


%%
%
% These are various ideas that I've tested while trying to get an iterative solver performance that is better then direct
% solver
%
% I have never been able to get the iterative solver to be even close to the direct one in terms of performance.
%
% The best option so far was use of gmres with ilu preconditioner combined with dissect. This was 'only' about 2 to 5 times slower than
% the direct solver. Also found the performance of the iterative solver to be highly problem dependent.
%
% Seems impossible to speed this up using parallel options as dissect not supported for distributed arrays (2023b)
% 
% If series of similar solves, then by only doing the ilu and dissect for the first solve, the following solves are faster than
% direct solve.
%
% Conclusions:
%
% ilu preconditioner using setup.type = "ilutp"; takes long time to produce, but convergence is good
%                    using 'nofill' is much faster, but convergence poor or no convergence
%
% gmres and bicgstab seem best, possibly bicgstab the better one
% 
% reordering using 'dissect' essential for memory and performance
%
% It appears that if ilu/ilutp pre-conditioner needs to be created, iterative solve is slower the direct. But if that
% preconditioner can be reused, e.g. in NR iteration, then iterative can be about 3 to 4 times faster using GPU. With CPU
% iterative approach is much slower than direct approach.
%
% Using single precision on GPU causes convergence issues and does not speed things up (most likely problem dependent).
%
% The surprising thing is that the ilu takes longer than a direct solve...?!
%
%%

NumWorkers=8 ;

ParPool = gcp('nocreate') ;

if isempty(ParPool)

    parpool('Processes',NumWorkers)

elseif (ParPool.NumWorkers~=NumWorkers)

    delete(gcp('nocreate'))
    parpool('Processes',NumWorkers)

end


% if isempty(ParPool)
% 
%     parpool('Threads',NumWorkers)
% 
% elseif (ParPool.NumWorkers~=NumWorkers)
% 
%     delete(gcp('nocreate'))
%     parpool('Threads',NumWorkers)
% 
% end



parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:genericError');
parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:SaveNotSupported');
warning('off','MATLAB:decomposition:genericError')
warning('off','MATLAB:decomposition:LoadNotSupported') 



%%

TestCase="-direct-" ;
% TestCase="-compare-" ;
TestCase="-best-"  ;


load("solveKApePIGTWGuvh250896.mat","A","B","CtrlVar","f","g","x0","y0")

CtrlVar.Parallel.Distribute=false;

switch TestCase


    case "-direct-"

        % Comparison with Direct solver

        
        tic
        [x,y,tolA,tolB]=ABfgPreEliminate(CtrlVar,A,B,f,g) ;
        toc
        x0=x;
        y0=y;


    case "-compare-"

        % Iterative method comparison


        Peq=[] ;
        Req=[] ;
        Ceq=[] ;
        L1=[];
        U1=[];
        L1eq=[];
        U1eq=[];

        

        load("solveKApePIGTWGuvh250896time0k19NRit2.mat","A","B","CtrlVar","f","g")
        CtrlVar.InfoLevelLinSolve=100;
        
        parpool('Threads')
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

%         From the method comparison it is concluded that ilutp+dissect+gmrs is the best approach.
%         Using equilibrate does not improve the convergence significantly, and takes lot of time (30 sec, for the typical example
%         used here).
% 
%         So fastest approach appears to be:
% 
%         1) Preconditioner based in ilutp with drop tolerance around 1e-6 or so
%         2) dissect  (a must to limit memory)
%         3) gmres with tol=1e-15 and maxit around 30 (Usually only five or so needed)
%
% 2025-05: It appears that iterative gmres is considerably faster on GPU than CPU.
% 
% It also seems that on GPU one can use a higher drop tolerance, which speeds up the ilu, at the cost of more iterations, and
% still be faster than the direct solver. 
%
%


        L=[]; U=[] ; P=[] ; x0=[] ; y0=[] ; perm=[] ; 


        % load solveKApePIGTWGuvh250896.mat ;
        %load("solveKApePIGTWGuvh250896time0k19NRit2.mat","A","B","CtrlVar","f","g","x0","y0")
        load("solveKApe843153_4083It1.mat","A","B","f","g","x0","y0")
        
        CtrlVar.ABfgPreEliminateIterative.Processor="-GPU-";
        CtrlVar.Parallel.Distribute=false;
        % A=distributed(A) ; B=distributed(B) ;  f=distributed(f) ; g=distributed(g) ; x0=distributed(x0) ; y0=distributed(y0) ;  % this does not work because dissect does not support distributed arrays
        % A=gpuArray(A) ; B=gpuArray(B) ;  f=gpuArray(f) ; g=gpuArray(g) ; x0=gpuArray(x0) ; y0=gpuArray(y0) ;  % this does not work because dissect does not support distributed arrays

        CtrlVar.InfoLevelLinSolve=100;  % 10 prints info, 100 creates figures as well

        x0=[] ; y0=[] ; xtilde0=[];
        fprintf("\n-----------------------------------------------------------------------------------------\n\n")

        fprintf(" Direct solve. \n")
        tDirect=tic;

        [xTest,yTest,dAtilde,tolA,tolB]=ABfgPreEliminate(CtrlVar,A,B,f,g) ;
        tDirect=toc(tDirect);
        fprintf(" Direct solve in %g \n",tDirect)
        fprintf("\t Tolerances for the direct solve are: %g \t %g \n",tolA,tolB)

        fprintf("\n-----------------------------------------------------------------------------------------\n\n")
        fprintf(" First iterative solve, involving creating pre-conditioner. \n")




        % solve iterative
        tstart1=tic ;
        [x,y,tolA,tolB,L,U,P,perm,xtilde]=ABfgPreEliminateIterative(CtrlVar,A,B,f,g,x0,y0,L,U,P,perm,xtilde0) ;
        tend1=toc(tstart1) ;



        fprintf(" Iterative solve in %g sec \t Direct solve in %g sec \n ",tend1,tDirect)
        fprintf("\t Tolerances for the iterative solve are: %g \t %g \n",tolA,tolB)


        fprintf("\n-----------------------------------------------------------------------------------------\n\n")
      
        % Maybe try loading a slightly difference matrix system, and use the previous pre-conditioner 
        %   load("solveKApePIGTWGuvh250896time0k19NRit3.mat","A","B","CtrlVar","f","g","x0","y0")
         fprintf(" loading another but similar system for the second iterative solve. \n") 
         load("solveKApe843153_4083It2.mat","A","B","f","g","x0","y0")
        % x0=x+1e-8*abs(x).*rand(length(x0),1) ; % for this slight modification of x0, a repeated solve is fast and does not require
        
        xtilde0=xtilde; % using previous solution as start
        % second iterative solve using preconditioner 
        fprintf(" Second iterative solve, now using previous pre-conditioner. \n")
        tstart2=tic ;
        [x,y,tolA,tolB,L,U,P,perm,xtilde]=ABfgPreEliminateIterative(CtrlVar,A,B,f,g,x0,y0,L,U,P,perm,xtilde0) ;
        tend2=toc(tstart2) ; 

       
        fprintf("\t Iterative solve in %g sec \t Direct solve in %g sec  \n ",tend2,tDirect)
        fprintf("\t Tolerances for the iterative solve are: %g \t %g \n",tolA,tolB)

     


    otherwise

        fprintf(" Which case? \n")

end



