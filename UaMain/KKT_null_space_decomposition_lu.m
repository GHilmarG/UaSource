


function [x, lambda] = KKT_null_space_decomposition_lu(H, g, A, b)


% KKT_NULL_SPACE_LU  Solves the KKT system using the null space method
%                    with sparse LU factorisation via the MATLAB
%                    decomposition object.
%
% Solves the KKT system:
%
%   | H   A' | | x      |   | g |
%   | A   0  | | lambda | = | b |
%
% Method:
%   1. Compute null space basis Z of A via sparse QR of A'
%   2. Compute particular solution xp satisfying A*xp = b
%   3. Form and factorise reduced system Hr = Z'*H*Z  ((n-m) x (n-m))
%   4. Solve Hr * u = Z'*(g - H*xp)  for u  (single right hand side)
%   5. Recover primal solution x = xp + Z*u
%   6. Recover multipliers lambda from A*A'*lambda = A*(g - H*x)
%
% Inputs:
%   H  - (n x n) sparse non-singular non-symmetric matrix
%   g  - (n x 1) right hand side vector
%   A  - (m x n) sparse full rank constraint matrix
%   b  - (m x 1) constraint right hand side
%
% Outputs:
%   x      - (n x 1) primal solution vector
%   lambda - (m x 1) dual solution vector (KKT multipliers)

% --- Input validation ---
[m, n] = size(A);


if isempty(A) || ( m==0 && n==0)   % In this case there is no A matrix and I just use the backslash operator. 
                  % This is here to cover this special case and make it possible to call this
                  % function even if there is not KKT system as such.

   x=H\g; lambda=[];

   return
end


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


is_AAT_Diag=isdiag(A*A');

if ~is_AAT_Diag
    dAAT = decomposition(A*A', 'chol');
    xp   = A' * (dAAT \ b);             % n x 1 particular solution
else

    xp   = A' * b;             % n x 1 particular solution
end

% --- Step 3: Form and factorise reduced system Hr = Z'*H*Z ---
% Hr is (n-m) x (n-m). Since H is non-symmetric, Hr is non-symmetric
% and LU factorisation is used instead of Cholesky.
Hr  = Z' * (H * Z);                  % (n-m) x (n-m)
dHr = decomposition(Hr, 'lu');

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
if ~is_AAT_Diag
    lambda = dAAT \ (A * (g - H * x));
else
    lambda = A * (g - H * x);
end


% % --- Diagnostic residuals (comment out in production) ---
% res1 = norm(H*x + A'*lambda - g)/(norm(x)+norm(lambda));
% res2 = norm(A*x - b)/norm(x);
% fprintf('Stationarity  ||H*x + A''*lambda - g||/(||x||+||lambda||)= %.2e\n', res1);
% fprintf('Feasibility   ||A*x - b||/||x||             = %.2e\n', res2);


end