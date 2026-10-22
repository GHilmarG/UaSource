


function [x, lambda] = KKT_null_space_decomposition_chol(H, g, A, b, CtrlVar, StructuralPattern)

%%
% KKT_NULL_SPACE_CHOL  Solves the KKT system using the null space method
%                      with sparse Cholesky factorisation via the MATLAB
%                      decomposition object.
%
% Solves the KKT system:
%
%   | H   A' | | x      |   | g |
%   | A   0  | | lambda | = | b |
%
% Method:
%   1. Compute null space basis Z of A via sparse QR of A'
%   2. Compute particular solution xp satisfying A*xp = b
%   3. Form and factorise reduced Hessian Hr = Z'*H*Z  ((n-m) x (n-m))
%   4. Solve Hr * u = Z'*(g - H*xp)  for u  (single right hand side)
%   5. Recover primal solution x = xp + Z*u
%   6. Recover multipliers lambda from A*A'*lambda = A*(g - H*x)
%   (16 Oct 2026) For nodal constraints (one nonzero in each row of A, all in different columns), steps 1-6 are replaced by a direct
%   reduction to the free unknowns, see below.
%
% Inputs:
%   H  - (n x n) sparse symmetric positive definite matrix
%   g  - (n x 1) right hand side vector
%   A  - (m x n) sparse full rank constraint matrix
%   b  - (m x 1) constraint right hand side
%
% Outputs:
%   x      - (n x 1) primal solution vector
%   lambda - (m x 1) dual solution vector (KKT multipliers)
%%



% --- Input validation ---
[m, n] = size(A);

% (16 Oct 2026) Optional fifth input CtrlVar: if CtrlVar.LinearSolver.Type="PARDISO", the system with H (if A is empty) and the
% system with the reduced Hessian Hr are solved with PARDISO (see UaPardisoSolve.m), the latter with a Cholesky factorisation.
UsePardiso = nargin>=5 && isfield(CtrlVar,"LinearSolver") && isfield(CtrlVar.LinearSolver,"Type") && CtrlVar.LinearSolver.Type=="PARDISO" ;
if nargin<6 ; StructuralPattern=[] ; end   % (17 Oct 2026) structural sparsity pattern of the mesh (eg MUA.uvAssemblyPattern), for PARDISO


if isempty(A) || ( m==0 && n==0)   % In this case there is no A matrix and I just use the backslash operator. 
    % This is here to cover this special case and make it possible to call this
    % function even if there is not KKT system as such.

    if UsePardiso
        x=UaPardisoSolve(CtrlVar,H,g,H,"symmetric",StructuralPattern); lambda=[];
    else
        x=H\g; lambda=[];
    end

    return
end


assert(n == size(H,1),  'H and A must have compatible dimensions');
assert(n == size(H,2),  'H must be square');
assert(m == length(b),  'A and b must have compatible dimensions');
assert(n == length(g),  'H and g must have compatible dimensions');
assert(m <= n,          'System must not be over-determined (m <= n)');

% --- (16 Oct 2026) Nodal constraints: direct reduction, without QR factorisation ---
% If each row k of A has exactly one nonzero, d_k in column c_k, and no column appears in more than one row, each constraint fixes one
% unknown: d_k x(c_k) = b_k. The null space of A is then spanned by the unit vectors of the free (unconstrained) unknowns f, ie Z=I(:,f),
% and the system can be reduced directly:
%
%     x(c) = b./d ,    H(f,f) x(f) = g(f) - H(f,c) x(c) ,    lambda = (g(c) - H(c,:) x)./d
%
% This is mathematically identical to the general null-space solution below (steps 1-6), but avoids the QR factorisation of A', the
% products Z'*H*Z and Z'*(g-H*xp), and the factorisation of A*A'. For all other forms of A (eg constraints that tie several unknowns
% together), the general method below is used.
[iRow, iCol, dVal] = find(A);
if numel(dVal)==m && numel(unique(iRow))==m && numel(unique(iCol))==m
    d=zeros(m,1) ; c=zeros(m,1) ;
    d(iRow)=dVal ; c(iRow)=iCol ;
    f=true(n,1) ; f(c)=false ;
    x=zeros(n,size(g,2)) ;
    x(c,:)=b./d+zeros(m,size(g,2)) ;
    if any(f)
        Hff=H(f,f) ;
        Hff=(Hff+Hff')/2 ;   % as Hr below: H is symmetric by construction, but may not be exactly so due to rounding errors
        rhs_f=g(f,:)-H(f,c)*x(c,:) ;
        if UsePardiso
            if isempty(StructuralPattern) ; SPf=[] ; else ; SPf=struct("Pattern",StructuralPattern,"Index",find(f)) ; end   % pattern of H(f,f)
            x(f,:)=UaPardisoSolve(CtrlVar,Hff,rhs_f,[],"spd",SPf) ;
        else
            x(f,:)=decomposition(Hff,'chol')\rhs_f ;
        end
    end
    lambda=(g(c,:)-H(c,:)*x)./d ;
    return
end

% --- Step 1: Compute null space basis Z of A via sparse QR of A' ---
% When called on a sparse matrix, qr automatically uses the spqr algorithm.
% Z is the last (n-m) columns of Q and spans the null space of A.
[Q, ~] = qr(A');
Z = Q(:, m+1:end);                   % n x (n-m) orthonormal null space basis

% --- Step 2: Compute particular solution xp satisfying A*xp = b ---
% Minimum norm solution via A'*(A*A')^{-1}*b.
% Factorise A*A' once and reuse in Step 6 for multiplier recovery.
dAAT = decomposition(A*A','chol');
%dAAT = decomposition(A*A');
xp   = A' * (dAAT \ b);             % n x 1 particular solution

% --- Step 3: Form and factorise reduced Hessian Hr = Z'*H*Z ---
% Hr is (n-m) x (n-m) and SPD if H is SPD on the null space of A
Hr  = Z' * (H * Z);                  % (n-m) x (n-m)
Hr=(Hr+Hr')/2 ; % symmetrize as the reduced Hessian must be symmetrical by construction, but may not always be so due to rounding errors. 
if ~UsePardiso
    dHr = decomposition(Hr, 'chol');
end

% --- Step 4: Solve reduced system Hr * u = Z'*(g - H*xp) ---
% Derived by substituting x = xp + Z*u into H*x + A'*lambda = g
% and premultiplying by Z' (the term Z'*A'*lambda vanishes since A*Z = 0)
rhs_r = Z' * (g - H * xp);
if UsePardiso
    u = UaPardisoSolve(CtrlVar,Hr,rhs_r,[],"spd") ;   % (16 Oct 2026) Cholesky factorisation with PARDISO
else
    u = dHr \ rhs_r;
end

% --- Step 5: Recover primal solution ---
x = xp + Z * u;

% --- Step 6: Recover multipliers lambda ---
% From A'*lambda = g - H*x, multiply both sides by A to get
% A*A'*lambda = A*(g - H*x), then reuse dAAT factorisation
lambda = dAAT \ (A * (g - H * x));


% 
% % --- Diagnostic residuals (comment out in production) ---
% res1 = norm(H*x + A'*lambda - g)/(norm(x)+norm(lambda)+eps);
% res2 = norm(A*x - b)/norm(x);
% fprintf('Stationarity  ||H*x + A''*lambda - g||/(||x||+||lambda||)= %g\n', res1);
% fprintf('            Feasibility   ||A*x - b||/||x||             = %g\n', res2);
% 



end















