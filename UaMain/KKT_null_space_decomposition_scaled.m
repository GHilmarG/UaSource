function [x, lambda] = KKT_null_space_decomposition_scaled(H, g, A, b)
% KKT_NULL_SPACE_SCALED  Solves the KKT system using the null space method
%                        with diagonal scaling and the MATLAB decomposition
%                        object.
%
% Solves the KKT system:
%
%   | H   A' | | x      |   | g |
%   | A   0  | | lambda | = | b |
%
% Diagonal scaling is applied to improve conditioning:
%   dx(i) = 1/sqrt(H(i,i))  -- scales diagonal of H to unity
%   dc(j) = 1/norm(A(j,:))  -- scales rows of A to unit 2-norm
%
% The scaled system is:
%   | Hs  As' | | xs      |   | gs |
%   | As  0   | | lambdas | = | bs |
%
% where:
%   Hs = Dx*H*Dx,  As = Dc*A*Dx,  gs = Dx*g,  bs = Dc*b
%
% After solving, the original variables are recovered as:
%   x = Dx * xs,  lambda = Dc * lambdas
%
% Method:
%   1. Compute scaling vectors dx and dc and form scaled system
%   2. Compute null space basis Z of As via sparse QR of As'
%   3. Compute particular solution xps satisfying As*xps = bs
%   4. Form and factorise reduced Hessian Hr = Z'*Hs*Z  ((n-m) x (n-m))
%   5. Solve Hr * u = Z'*(gs - Hs*xps)  for u  (single right hand side)
%   6. Recover scaled primal solution xs = xps + Z*u
%   7. Recover scaled multipliers lambdas from As*As'*lambdas = As*(gs-Hs*xs)
%   8. Unscale to recover x and lambda
%
% Inputs:
%   H  - (n x n) sparse non-singular matrix
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

% --- Step 1: Compute scaling vectors and form scaled system ---
% dx: scale so diagonal of H is unity
% dc: scale so rows of A have unit 2-norm
%dx = 1 ./ sqrt(abs(diag(H)));        % n x 1
dx = 1./sqrt(sum(H.^2,2));
dc = 1 ./ sqrt(sum(A.^2, 2));        % m x 1

% Form sparse diagonal scaling matrices
Dx = spdiags(dx, 0, n, n);
Dc = spdiags(dc, 0, m, m);

% Scaled system matrices and right hand sides
Hs = Dx * H * Dx;                    % n x n
As = Dc * A * Dx;                    % m x n
gs = Dx * g;                          % n x 1
bs = Dc * b;                          % m x 1

% --- Step 2: Compute null space basis Z of As via sparse QR of As' ---
% When called on a sparse matrix, qr automatically uses the spqr algorithm.
% Z is the last (n-m) columns of Q and spans the null space of As.
[Q, ~] = qr(As');
Z = Q(:, m+1:end);                   % n x (n-m) orthonormal null space basis

% --- Step 3: Compute particular solution xps satisfying As*xps = bs ---
% Minimum norm solution via As'*(As*As')^{-1}*bs.
% Factorise As*As' once and reuse in Step 7 for multiplier recovery.
dAAsT = decomposition(As*As', 'auto');
xps   = As' * (dAAsT \ bs);         % n x 1 particular solution

% --- Step 4: Form and factorise reduced Hessian Hr = Z'*Hs*Z ---
% Hr is (n-m) x (n-m)
Hr  = Z' * (Hs * Z);                 % (n-m) x (n-m)
dHr = decomposition(Hr, 'auto');

% --- Step 5: Solve reduced system Hr * u = Z'*(gs - Hs*xps) ---
% Only a single right hand side regardless of m
rhs_r = Z' * (gs - Hs * xps);
u     = dHr \ rhs_r;

% --- Step 6: Recover scaled primal solution ---
xs = xps + Z * u;

% --- Step 7: Recover scaled multipliers ---
% From As'*lambdas = gs - Hs*xs, multiply both sides by As to get
% As*As'*lambdas = As*(gs - Hs*xs), then reuse dAAsT factorisation
lambdas = dAAsT \ (As * (gs - Hs * xs));

% --- Step 8: Unscale ---
x      = Dx * xs;
lambda = Dc * lambdas;

% --- Diagnostic residuals (comment out in production) ---
res1 = norm(H*x + A'*lambda - g)/norm(x);
res2 = norm(A*x - b)/norm(x);
fprintf('Stationarity  ||H*x + A''*lambda - g||/||x|| = %.2e\n', res1);
fprintf('Feasibility   ||A*x - b||/||x||              = %.2e\n', res2);

end