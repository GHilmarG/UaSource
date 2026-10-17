function [x,y,dAtilde,tolA,tolB]=ABfgPreEliminate(CtrlVar,A,B,f,g,dAtilde,StructuralPattern)


narginchk(5,7)

if nargin< 6
    dAtilde=[];
end
if nargin< 7
    StructuralPattern=[];   % (17 Oct 2026) structural sparsity pattern of the mesh, for PARDISO (see UaPardisoSolve.m)
end

% (16 Oct 2026) PARDISO (Intel oneMKL) instead of decomposition/mldivide, if CtrlVar.LinearSolver.Type="PARDISO" (see UaPardisoSolve.m)
UsePardiso = isfield(CtrlVar,"LinearSolver") && isfield(CtrlVar.LinearSolver,"Type") && CtrlVar.LinearSolver.Type=="PARDISO" ;

[nA,mA]=size(A);
[nB,mB]=size(B);
[nf,mf]=size(f);


if isempty(B) && isempty(g) && ~isempty(A) && ~isempty(f) && mA==nf

    % Possibly not needed, but check if this is not just a very simple case of B=g=[]

    if UsePardiso
        x=UaPardisoSolve(CtrlVar,A,f,A,"nonsymmetric",StructuralPattern) ;
    else
        x=A\f;
    end
    y=NaN;

    if nargout>3
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
        if isempty(dAtilde)      % (15 Oct 2026) Atilde is only needed for the decomposition, which may be given
            Atilde=Q*A+ BtB ;
        else
            Atilde=[] ;
        end
        btilde=(Q*f+B'*g) ; % when solving for multiple right-hand sides, this will use automatic implicit expansion to expand B'*g to match Q*f

  

        if CtrlVar.Parallel.isTest && ~isempty(Atilde)

            % seq:

            % distribute
            tSeq=tic ;
            dAtildeSeq=decomposition(Atilde);  % this might not be needed if this has been done already, but for speed comparison this is done here each time
            xSeq=dAtildeSeq\btilde;
            tSeq=toc(tSeq) ; 


            % distribute
            tDistributed=tic ;
            AtildeDist=distributed(Atilde);
            btildeDist=distributed(btilde);
            dAtildeDist=decomposition(AtildeDist);  % this might not be needed if this has been done already, but for speed comparison this is done here each time
            xDist=dAtildeDist\btildeDist;
            tDistributed=toc(tDistributed) ; 
        
            [nAtilde,mAtilde]=size(Atilde);
            fprintf('\n ----------------------------- Info on distributed solve performance : \n')
            fprintf('%i x %i :  tSeq=%f \t tDistributed=%f \t speedup=%f \n',nAtilde,mAtilde,tSeq,tDistributed,tSeq/tDistributed) ;
            fprintf(' norm(xSeq-xDist)/norm(xSeq)=%g    \n',full(norm(xSeq-xDist)/norm(xSeq)))
            fprintf(' ----------------------------- \n')


        end


        % https://uk.mathworks.com/help/parallel-computing/benchmarking-a-b.html

        
        if CtrlVar.Parallel.Distribute && ~isempty(Atilde)
            if ~isdistributed(Atilde)
                Atilde=distributed(Atilde);
            end
            if ~isdistributed(btilde)
                btilde=distributed(btilde);
            end
        end

        %  decomposition is about the same, and as expected this only speeds things up if several solves with the same matrix
        %  are needed.
        %
        tDecomposition=tic;
        if UsePardiso && ~isdistributed(Atilde)
            % (16 Oct 2026) factorisation and solve with PARDISO, see UaPardisoSolve.m. The analysis (reordering and symbolic
            % factorisation) is kept between calls, and only repeated if the sparsity pattern of Atilde changes.
            if ~isempty(dAtilde)
                error("ABfgPreEliminate:PARDISO","Reusing a decomposition (modified Newton-Raphson) is not supported with CtrlVar.LinearSolver.Type=""PARDISO"".")
            end
            x=UaPardisoSolve(CtrlVar,Atilde,btilde,A,"nonsymmetric",StructuralPattern) ;   % the patterns of A and of the mesh contain all rows of Atilde
            tDecomposition=toc(tDecomposition) ;
            tSolve=0 ;
        else
            if isempty(dAtilde)
                dAtilde=decomposition(Atilde);
            end
            tDecomposition=toc(tDecomposition) ;

            tSolve=tic;
            x=dAtilde\btilde;
            tSolve=toc(tSolve);
        end

        % tSolve=tic;
        % x=Atilde\btilde;
        % tSolve=toc(tSolve);

        % [tDecomposition tSolve]

        if isdistributed(x)
            x=gather(x) ;
        end
        % toc


     







        
        y=B*(f-A*x);

        % Now the solution of the scaled system has been found.

        y=y/factor;

        if nargout>3
            % check if within tolerances

            A=A/factor ;
            f=f/factor ;
            tolA=norm(A*x+B'*y-f)/norm(f);
            tolB=norm(B*x-g);

            if tolA>1e-6 || tolB>1e-6

                fprintf('ABfgPreEliminate: Error seems too large or \t \t \t %g \t %g \n ',norm(A*x+B'*y-f)/norm(f),norm(B*x-g))

            end
        end

        y=Scale*y; % and now scale y in case B and g were scaled above.

    else


        error('ABfgPreEliminate:B','B*B^T not diagonal')
    end
    % % --- Diagnostic residuals (comment out in production) ---
    % if nargout<=3
    %     A=A/factor ;
    %     f=f/factor ;
    % end
    % res1 = norm(A*x + B'*y - f)/(norm(x)+norm(y));
    % res2 = norm(B*x - g)/norm(x);
    % fprintf('Stationarity  ||H*x + A''*lambda - g||/(||x||+||lambda||)= %.2e\n', res1);
    % fprintf('Feasibility   ||A*x - b||/||x||             = %.2e\n', res2);
end

end

