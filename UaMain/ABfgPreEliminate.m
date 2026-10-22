function [x,y,dAtilde,tolA,tolB]=ABfgPreEliminate(CtrlVar,A,B,f,g,dAtilde,StructuralPattern)

%%
%
% I developed this approach independently, and then later learned that this exact method is Eq (29) in
%
% F. Aghili, "A unified approach for inverse and direct dynamics of constrained multibody systems based 
% on linear projection operator: applications to control and simulation",
% IEEE Trans. on Robotics 21(5), 834–849 (2005), doi:10.1109/TRO.2005.851380.
%
% An almost identical approach is also found in:
%
% M. Ainsworth, "Essential boundary conditions and multi-point constraints in finite element analysis", 
% Computer Methods in Applied Mechanics and Engineering, 190 (48). pp. 6323-6339.
%
%
%
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
[nf,~]=size(f);


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
    
        if UsePardiso && ~isdistributed(Atilde)
            % (16 Oct 2026) factorisation and solve with PARDISO, see UaPardisoSolve.m. The analysis (reordering and symbolic
            % factorisation) is kept between calls, and only repeated if the sparsity pattern of Atilde changes.
            if ~isempty(dAtilde)
                error("ABfgPreEliminate:PARDISO","Reusing a decomposition (modified Newton-Raphson) is not supported with CtrlVar.LinearSolver.Type=""PARDISO"".")
            end
            x=UaPardisoSolve(CtrlVar,Atilde,btilde,A,"nonsymmetric",StructuralPattern) ;   % the patterns of A and of the mesh contain all rows of Atilde
        
        else
            if isempty(dAtilde)
                dAtilde=decomposition(Atilde);
            end
     

        
            x=dAtilde\btilde;
         
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

