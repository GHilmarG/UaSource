

function [x, lambda] = KKT_range_space_symmetric(H, g, A, b)


% KKT_RANGE_SPACE_CHOL  Solves the KKT system using the range space method
%                       with sparse Cholesky factorisation.
%
% Solves the KKT system:
%
%   | H   A' | | x      |   | g |
%   | A   0  | | lambda | = | b |
%
% Method:
%   1. Factorise H using sparse Cholesky with vector permutation
%   2. Form Schur complement Sc = A * H^{-1} * A'  (m x m)
%   3. Solve Sc * lambda = A * H^{-1} * g - b  for lambda
%   4. Recover x = H^{-1} * (g - A' * lambda)
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

% --- Step 1: Sparse Cholesky factorisation of H ---
% Factorisation satisfies: H(p,p) = R'*R
% Solution to H * x = rhs recovered as: x(p) = R \ (R' \ rhs(p))
[R, flag, p] = chol(H, 'vector');
assert(flag == 0, 'H is not positive definite');

% --- Step 2: Form Schur complement Sc = A * H^{-1} * A' ---
HinvAt        = sparse(n, m);
HinvAt(p,:)   = R \ (R' \ A(:,p)');
Sc            = A * HinvAt;          % m x m Schur complement

% --- Step 3: Solve for multipliers lambda ---
% A * H^{-1} * A' * lambda = A * H^{-1} * g - b
Hinvg    = zeros(n, 1);
Hinvg(p) = R \ (R' \ g(p));
lambda   = Sc \ (A * Hinvg - b);

% --- Step 4: Recover primal solution x ---
% x = H^{-1} * (g - A' * lambda)
rhs  = g - A' * lambda;
x    = zeros(n, 1);
x(p) = R \ (R' \ rhs(p));

end