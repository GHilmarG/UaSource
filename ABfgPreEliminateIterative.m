





function [x,y,tolA,tolB,L,U,P,perm,xtilde]=ABfgPreEliminateIterative(CtrlVar,A,B,f,g,x0,y0,L,U,P,perm,xtilde0)

%%
%
%
%
%

nargoutchk(9,9)
narginchk(12,12)

persistent iFigure

if isempty(iFigure)
    iFigure=100;
end

if nargin < 8
    L=[];
    U=[];

end

if nargin<6
    x0=[];
    y0=[];
end


[nA,mA]=size(A);
[nB,mB]=size(B);
[nf,mf]=size(f);

if isempty(x0)
    x0=f*0;
end

if isempty(xtilde0)
    xtilde=f*0;
end





% ilutp results in much smaller residuals at fewer iterations, but ilu takes long time

setup.type = "nofill"; setup.milu = "off"; 
tol=1e-13 ; maxit=5; restart=50;   % quick for testing purposes

% takes forever
% setup.type = "crout"; setup.milu = "off"; setup.droptol = 1e-6 ; 
% tol=1e-13 ; maxit=5; restart=5;   % quick for testing purposes

setup.type = "ilutp"; setup.milu = "off"; setup.droptol = 1e-6;    setup.udiag=0 ;  
tol=1e-13 ; maxit=2; restart=50;   % quick for testing purposes


if isempty(B) && isempty(g) && ~isempty(A) && ~isempty(f) && mA==nf

    % Possibly not needed, but check if this is not just a very simple case of B=g=[]

    x=A\f;
    y=NaN;

    if nargout>2
        tolB=NaN;
        tolA=norm(A*x-f)/norm(f);
    end

else

    BBT=B*B';

    if isdiag(BBT)  % the method assumes that B B' is diagonal


        tCPUtotal=tic;

        % It also assumes that B B' is a unity matrix, but if not then simple scaling can be used
        % to ensure that this is the case.
        % To make this a bit more general, I here check if B B' is indeed unity, and
        % if not I do the requried scaling.
        tolerance=eps*1000;
        isBBTunity=all(abs(diag(BBT) - 1) < tolerance) ;

        if ~isBBTunity
            [B,g,~,Scale]=ScaleL(CtrlVar,B,g) ;
        else
            Scale=1;
        end

        % For numerical reasons a further simple scaling of A is done to bring the
        % sizes of the elements of A in line with those of B.
        factor=1./((mean(full(abs(diag(A))))));  % mean does not work for gpuArray !!!

        if ~isfinite(factor)  % just in case all elements along the diagonal happen to be equal to zero
            factor=1;
        end

        BtB=B'*B;
        A=factor*A ; f=factor*f ;  % this leaves x unaffected but y is scaled
        Q=speye(nA,nA)-BtB ;
        Atilde=Q*A+ BtB ;
        btilde=(Q*f+B'*g) ;


        %% ilu for both equilibrated and not


        tdissectAtilde=tic;
        if isempty(perm)  % If the matrix has the same sparsity structure, then I don't need to do the permutation again
            perm=dissect(Atilde);  % does not work for distributed or gpuarrays

            % The nested dissection algorithm produces high quality reordering and performs particularly well with finite element
            % matrices compared to other reordering techniques. For more information about the nested dissection ordering
            % algorithm, see https://uk.mathworks.com/help/matlab/math/sparse-matrix-operations.html

        end
        % pAtilde=symrcm(Atilde);      % ilu appears to take longer as compared to using dissect
        tdissectAtilde=toc(tdissectAtilde);

        iperm(perm)=1:length(perm);   % inverse of the permutation vector
        Atilde=Atilde(perm,perm) ;
        btilde=btilde(perm);
        x0=x0(perm) ; 
        % If I return xtilde and the re-use the same ordering (i.e. perm) then I can take the previous xtilde as a starting point,
        % and I don't need to re-order that initial guess.

        % can't use setup.type = "ilutp" with distributed arrays
        % but even for setup.type = "ilutp"; setup.milu = "off"; setup.droptol = 1e-6;    setup.udiag=0 ;
        % using distributed arrays is slower...
        tluinc=tic;
        if isempty(L)
            [L,U,P] = ilu(Atilde,setup);
        end
        tluinc=toc(tluinc);

        if setup.type == "ilutp"
            Atilde=P*Atilde ; btilde=P*btilde;
        end
       
        %fprintf("norm(Atilde-L*U)/norm(Atilde)=%g \n",norm(Atilde-L*U,"fro")/norm(Atilde,"fro"))

   %     [sol,flag,relres,iter,resvec]=bicgstabl(Atilde,btilde,tol,maxit,L,U,xtilde0);


        tgmres=tic;

        [xtilde,flag,relres,iter,resvec]=gmres(Atilde,btilde,restart,tol,maxit,L,U,xtilde0);

        x=xtilde(iperm) ;
        tgmres=toc(tgmres);

        % Important to replace zeros on diagonal, but does not converge well, much worse that using ilu
        % k = 3;
        % M = tril(triu(Atilde,-k),k);
        % I=find(abs(diag(M))==0 );   % must make sure to replace 0 on the diagonal with 1
        % [n,m]=size(Atilde) ;
        % M=M+sparse(I,I,1,n,m); 
        % [xtilde,flag,relres,iter,resvec]=gmres(Atilde,btilde,restart,tol,maxit,M);
        % 



        % % Not particularly fast, and less good convergence, possibly because not possible to provide L and U separately 
        % fprintf(" GPU \n ")
        % tic
        % AtildeGPU=gpuArray(Atilde) ; M=L*U ; MGPU=gpuArray(M) ; btildeGPU=gpuArray(btilde) ;
        % toc
        % 
        % tic
        % [xGPU,flag,relresGPU,iter,resvecGPU]=gmres(AtildeGPU,btildeGPU,restart,tol,maxit,MGPU,[],xtilde0);
        % toc
        % % 



        % For some reason, very slow
        % tic
        % AtildeDist=distributed(Atilde) ; Ldist=distributed(L) ; Udist=distributed(U) ; x0dist=distributed(x0) ; 
        % [x,flag,relres,iter,resvec]=gmres(AtildeDist,btilde,restart,tol,maxit,Ldist,Udist,x0dist);
        % toc


        % fprintf(" spmd \n ")
        % For some reason, very slow
        % tic
        % % AtildeDist=distributed(Atilde) ; Ldist=distributed(L) ; Udist=distributed(U) ; x0dist=distributed(x0) ;
        % spmd
        %     AtildeDist=codistributed(Atilde); 
        %     [x,flag,relres,iter,resvec]=gmres(AtildeDist,btilde,restart,tol,maxit,L,U,x0);
        % end
        % toc


        tCPUtotal=toc(tCPUtotal) ;

        if CtrlVar.InfoLevelLinSolve>=10

            fprintf("                   dissect Atilde %f sec\n",tdissectAtilde)
            fprintf("                              ilu %f sec\n",tluinc)
            fprintf("                            gmres %f sec\n",tgmres)
            fprintf("total time for iterative solution %f sec\n",tCPUtotal)

            figure(iFigure) ; iFigure=iFigure+100;
           
            fprintf('\n flag=%-i, iter=%-i %-i, relres=%-g \n ',flag,iter(1),iter(2),relres)

            nnn=numel(resvec);
            semilogy(0:nnn-1,resvec,'-o',LineWidth=2)
            xlabel('Iteration Number',Interpreter='latex')
            ylabel('Relative Residual',Interpreter='latex')


            % fig = gcf; exportgraphics(fig,'IterativeSolveExample.pdf')
        end


        fprintf("\n\n")


        % x=Atilde\btilde;
        y=B*(f-A*x);

        % Now the solution of the scaled system has been found.

        y=y/factor;


        if nargout>2
            % check if within tolerances

            A=A/factor ;
            f=f/factor ;

            tolA=norm(A*x+B'*y-f)/norm(f);
            tolB=norm(B*x-g);

            %                if tolA>1e-6 || tolB>1e-6

            fprintf('\t residuals \t %g \t %g \n ',tolA,tolB)
            %
        end


        y=Scale*y; % and now scale y in case B and g were scaled above.



    else


        error('ABfgPreEliminate:B','B*B^T not diagonal')
    end

end




