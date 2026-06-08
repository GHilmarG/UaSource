


function [x, lambda] = KKT_null_space_decomposition_chol(H, g, A, b)
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

% --- Input validation ---
[m, n] = size(A);
assert(n == size(H,1),  'H and A must have compatible dimensions');
assert(n == size(H,2),  'H must be square');
assert(m == length(b),  'A and b must have compatible dimensions');
assert(n == length(g),  'H and g must have compatible dimensions');
assert(m <= n,          'System must not be over-determined (m <= n)');

% --- Step 1: Compute null space basis Z of A via sparse QR of A' ---
% When called on a sparse matrix, qr automatically uses the spqr algorithm.
% Z is the last (n-m) columns of Q and spans the null space of A.
[Q, ~] = qr(A');
Z = Q(:, m+1:end);                   % n x (n-m) orthonormal null space basis

% --- Step 2: Compute particular solution xp satisfying A*xp = b ---
% Minimum norm solution via A'*(A*A')^{-1}*b.
% Factorise A*A' once and reuse in Step 6 for multiplier recovery.
%dAAT = decomposition(A*A','chol');
dAAT = decomposition(A*A');
xp   = A' * (dAAT \ b);             % n x 1 particular solution

% --- Step 3: Form and factorise reduced Hessian Hr = Z'*H*Z ---
% Hr is (n-m) x (n-m) and SPD if H is SPD on the null space of A
Hr  = Z' * (H * Z);                  % (n-m) x (n-m)
Hr=(Hr+Hr')/2 ; % symmetrize as the reduced Hessian must be symmetrical by construction, but may not always be so due to rounding errors. 
dHr = decomposition(Hr, 'chol');

% --- Step 4: Solve reduced system Hr * u = Z'*(g - H*xp) ---
% Derived by substituting x = xp + Z*u into H*x + A'*lambda = g
% and premultiplying by Z' (the term Z'*A'*lambda vanishes since A*Z = 0)
rhs_r = Z' * (g - H * xp);
u     = dHr \ rhs_r;

% --- Step 5: Recover primal solution ---
x = xp + Z * u;

% --- Step 6: Recover multipliers lambda ---
% From A'*lambda = g - H*x, multiply both sides by A to get
% A*A'*lambda = A*(g - H*x), then reuse dAAT factorisation
lambda = dAAT \ (A * (g - H * x));



% --- Diagnostic residuals (comment out in production) ---
res1 = norm(H*x + A'*lambda - g)/(norm(x)+norm(lambda)+eps);
res2 = norm(A*x - b)/norm(x);
fprintf('Stationarity  ||H*x + A''*lambda - g||/(||x||+||lambda||)= %g\n', res1);
fprintf('            Feasibility   ||A*x - b||/||x||             = %g\n', res2);




end















