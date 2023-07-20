function [x,y,tolA,tolB,L,U,perm]=ABfgPreEliminateIterative(CtrlVar,A,B,f,g,x0,y0,L,U,perm)

%%
%
%
%
%

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


setup.type = 'nofill'; setup.milu = 'off'; 
setup.type = "ilutp"; setup.milu = "off"; setup.droptol = 1e-6;    setup.udiag=0 ;  isReorder=true; % must be used with re-ordering

tol=1e-13 ; maxit=2; restart=10;   % quick for testing purposes


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
        if isempty(perm)
            perm=dissect(Atilde);  % does not work for disstributed or gpuarrays 
        end
        % pAtilde=symrcm(Atilde);      % ilu appears to take longer as compared to using dissect
        tdissectAtilde=toc(tdissectAtilde);

        iperm(perm)=1:length(perm);   % inverse of the permutation vector
        Atilde=Atilde(perm,perm) ;
        btilde=btilde(perm);
        x0=x0(perm) ; 

        tluinc=tic;
        if isempty(L)
            [L,U] = ilu(Atilde,setup);
        end
        tluinc=toc(tluinc);


        %[sol,flag,relres,iter,resvec]=bicgstabl(AA,ff,tol,maxit,L1,U1,xx0);


        tgmres=tic;



        tic

        [x,flag,relres,iter,resvec]=gmres(Atilde,btilde,restart,tol,maxit,L,U,x0);
        

        toc

%         tic
%         AtildeGPU=gpuArray(Atilde) ; M=L*U ; MGPU=gpuArray(M) ; btildeGPU=gpuArray(btilde) ;
%         toc
% 
%         tic
%         [xGPU,flag,relresGPU,iter,resvecGPU]=gmres(AtildeGPU,btildeGPU,restart,tol,maxit,MGPU,[],x0);
%         toc

        x=x(iperm) ;
        tgmres=toc(tgmres);

        if CtrlVar.InfoLevelLinSolve>=10

            fprintf("  dissect Atilde %f sec\n",tdissectAtilde)
            fprintf("             ilu %f sec\n",tluinc)
            fprintf("          gmres %f sec\n",tgmres)
            fprintf("total time for nonEq with  is about %f",tdissectAtilde+tluinc+tgmres)

            figure
            fprintf("\n\n")
            fprintf(' flag=%-i, iter=%-i %-i, relres=%-g \n ',flag,iter(1),iter(2),relres)

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




