




function [x,y,tolA,tolB,Peq,Req,Ceq,L1,U1,L1eq,U1eq]=ABfgPreEliminateIterativeMethodComparision(CtrlVar,A,B,f,g,x0,y0,Peq,Req,Ceq,L1,U1,L1eq,U1eq)

%% This is a test run
%
% The purpose is to run various iterative options for comparision
%
%
%

if nargin < 8
    Peq=[];
    Req=[];
    Ceq=[];
    L1=[];
    U1=[];
    L1eq=[];
    U1eq=[];
end

if nargin<6
    x0=[]; 
    y0=[];
end


[nA,mA]=size(A);
[nB,mB]=size(B);
[nf,mf]=size(f);


if ~isempty(x0) && ~isempty(y0)

    tolA=norm(A*x0+B'*y0-f)/norm(f);
    tolB=norm(B*x0-g);
    fprintf(" Initial normed residuals are: \n")
    fprintf("   norm(A*x+B'*y-f)/norm(f) = %g \n",tolA) 
    fprintf("                  norm(B*x-g)= %g \n",tolB)

end


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
        % if not I do the required scaling.
        tolerance=eps*1000;
        isBBTunity=all(abs(diag(BBT) - 1) < tolerance) ;

        if ~isBBTunity
            [B,g,~,Scale]=ScaleL(CtrlVar,B,g) ;
        else
            Scale=1;
        end

        % For numerical reasons a further simple scaling of A is done to bring the
        % sizes of the elements of A in line with those of B.
        factor=1./(full(mean(abs(diag(A)))));

        if ~isfinite(factor)  % just in case all elements along the diagonal happen to be equal to zero
            factor=1;
        end

        BtB=B'*B;

        A=factor*A ; f=factor*f ;  % this leaves x unaffected but y is scaled



        Q=speye(nA,nA)-BtB ;
        Atilde=Q*A+ BtB ;
        btilde=(Q*f+B'*g) ;


        teq=tic;
        % c1 = condest(Atilde)
        %% Equilibriate  (this seems to take surprisingly long time, better to save locally when testing using same matrix)
        if isempty(Peq)
            [Peq,Req,Ceq] = equilibrate(Atilde);
        end
        AtildeEquilibrated = Req*Peq*Atilde*Ceq;
        btildeEquilibrated = Req*Peq*btilde;


        teq=toc(teq) ;


        %% ilu for both equilibrated and not


        setup.type = 'nofill'; setup.milu = 'off';  isReorder=false;
        setup.type = "ilutp"; setup.milu = "off"; setup.droptol = 1e-8;    setup.udiag=0 ;  isReorder=true; % must be used with re-ordering

        

        if isReorder

            tdissectAtilde=tic;
            pAtilde=dissect(Atilde);
            % pAtilde=symrcm(Atilde);      % ilu appears to take longer as compared to using dissect
            tdissectAtilde=toc(tdissectAtilde);

            tdissectBeq=tic;
            pAtildeEquilibrated=dissect(AtildeEquilibrated);
            % pBeq=symrcm(Beq);
            tdissectBeq=toc(tdissectBeq);

            invpAtilde(pAtilde)=1:length(pAtilde);   % inverse of the permutation vector
            invpBeq(pAtildeEquilibrated)=1:length(pAtildeEquilibrated);            % inverse of the permutation vector

            Atilde=Atilde(pAtilde,pAtilde) ;
            btilde=btilde(pAtilde);

            AtildeEquilibrated=AtildeEquilibrated(pAtildeEquilibrated,pAtildeEquilibrated) ;
            btildeEquilibrated=btildeEquilibrated(pAtildeEquilibrated);
        else

            tdissectAtilde=0;
            tdissectBeq=0; 

        end

        tluincEq=tic;
        if isempty(L1eq)
            [L1eq,U1eq] = ilu(AtildeEquilibrated,setup);
        end
        tluincEq=toc(tluincEq);

        tluinc=tic;
        if isempty(L1)
            [L1,U1] = ilu(Atilde,setup);
        end
        tluinc=toc(tluinc);



        tol=1e-15 ; maxit=1; restart=300;
        tol=1e-15 ; restart=10;  maxit=floor(15/restart) ; % quick for testing purposes



        %[sol,flag,relres,iter,resvec]=bicgstabl(AA,ff,tol,maxit,L1,U1,xx0);

        if isempty(x0)
            x0=f*0;
        end

        % gmres without preconditioners
        tgmres=tic;
        [x,flag,relres,iter,resvec]=gmres(Atilde,btilde,restart,tol,maxit,[],[],x0);

        if isReorder
            x=x(invpAtilde) ;
        end
        tgmres=toc(tgmres);


        % gmres with il0 preconditioners
        tgmresPre=tic;
        [xPre,flagPre,relresPre,iterPre,resvecPre]=gmres(Atilde,btilde,restart,tol,maxit,L1,U1,x0);
        if isReorder
            xPre=xPre(invpAtilde) ;
        end
        tgmresPre=toc(tgmresPre);

         % gmres, equilibrated matrix but no preconditioner
        xEq0=f*0;
        tgmresEq=tic;
        [xeq,flagEq,relresEq,iterEq,resvecEq]=gmres(AtildeEquilibrated,btildeEquilibrated,restart,tol,maxit,[],[],xEq0);
        if isReorder
            xeq=xeq(invpBeq) ;
        end
        xeq=Ceq*xeq;
        tgmresEq=toc(tgmresEq);

        % gmres, equilibrated matrix with preconditioners
        tgmresEqPre=tic;
        [xeqPre,flagEqPre,relresEqPre,iterEqPre,resvecEqPre]=gmres(AtildeEquilibrated,btildeEquilibrated,restart,tol,maxit,L1eq,U1eq,xEq0);
        if isReorder
            xeqPre=xeqPre(invpBeq) ;
        end
        xeqPre=Ceq*xeqPre;
        tgmresEqPre=toc(tgmresEqPre);


        sol{1}=x;      text(1)="x     " ;
        sol{2}=xPre;   text(2)="xPre  ";
        sol{3}=xeq;    text(3)="eq    ";
        sol{4}=xeqPre; text(4)="xeqPre" ;


        

        if CtrlVar.InfoLevelLinSolve>=1
            

            fprintf("      equlibrate %f sec\n",teq)


            fprintf("  dissect Atilde %f sec\n",tdissectAtilde)
            fprintf("    dissect ABEq %f sec\n",tdissectBeq)

            fprintf("             ilu %f sec\n",tluinc)
            fprintf("           iluEq %f sec\n",tluincEq)

            fprintf("          gmres %f sec\n",tgmres)
            fprintf("       gmresPre %f sec\n",tgmresPre)
            fprintf("        gmresEq %f sec\n",tgmresEq)
            fprintf("     gmresEqPre %f sec\n",tgmresEqPre)

            fprintf("\n total time for nonEq with pre is about %f",tdissectAtilde+tluinc+tgmresPre)



        end

        if CtrlVar.InfoLevelLinSolve>=10

            figure
            fprintf("\n\n")
            fprintf(' flag=%-i, iter=%-i %-i, relres=%-g \n ',flag,iter(1),iter(2),relres)
            fprintf(' flag=%-i, iter=%-i %-i, relresPre=%-g \n ',flagPre,iterPre(1),iterPre(2),relresPre)
            fprintf(' flag=%-i, iter=%-i %-i, relresEq=%-g \n ',flagEq,iterEq(1),iterEq(2),relresEq)
            fprintf(' flag=%-i, iter=%-i %-i, relresEqPre=%-g \n ',flagEqPre,iterEqPre(1),iterEqPre(2),relresEqPre)

            nnn=numel(resvec);
            semilogy((0:nnn-1)/2,resvec,'-',LineWidth=2)
            hold on
            nnn=numel(resvecPre);
            semilogy((0:nnn-1)/2,resvecPre,'-',LineWidth=2)
            nnn=numel(resvecEq);
            semilogy((0:nnn-1)/2,resvecEq,'-',LineWidth=2)
            nnn=numel(resvecEqPre);
            semilogy((0:nnn-1)/2,resvecEqPre,'-',LineWidth=2)
            xlabel('Iteration Number',Interpreter='latex')
            ylabel('Relative Residual',Interpreter='latex')
            legend("Original system","With ilu(0) preconditioning","Equlibrated","Equlibrated and ilu(0) preconditioning",interpreter="latex")

            % fig = gcf; exportgraphics(fig,'IterativeSolveExample.pdf')
        end


        fprintf("\n\n")
        TolAmin=inf;
        isScaled=true; 
        for I=1:4   % now plug the four solutions into the system and estimate norm of residuals


            xx=sol{I} ;

            % x=Atilde\btilde;

            if isScaled

                A=A/factor ;
                f=f/factor ;

                B=B/Scale;
                g=g/Scale;
                isScaled=false;

            end


            % yy=B*(f-A*xx);  %   B' yy= (f-A*xx); but since B B'= 1 this is correct, but then must use the scaled B here

            yy= ( B*B') \ (B*(f-A*xx)) ;

            tolA=norm(A*xx+B'*yy-f)/norm(f);
            tolB=norm(B*xx-g);

            %                if tolA>1e-6 || tolB>1e-6

            fprintf('\t %s \t %g \t %g \n ',text(I),tolA,tolB)
            %

            % Select the solution with the lowers tolerance as the one returned

   %         yy=Scale*yy; % and now scale y in case B and g were scaled above.

            if tolA < TolAmin
                x=xx;
                y=yy;
                TolAmin=tolA;
            end


        end

        



    else


        error('ABfgPreEliminate:B','B*B^T not diagonal')
    end

end




