


function [x, lambda] = KKT_range_space_unsymmetric(H, g, A, b)

%%
% KKT_RANGE_SPACE_LU  Solves the KKT system using the range space method
%                     with sparse LU factorisation.
%
% Solves the KKT system:
%
%   | H   A' | | x      |   | g |
%   | A   0  | | lambda | = | b |
%
% Method:
%   1. Factorise H using sparse LU with vector permutations
%   2. Form Schur complement Sc = A * H^{-1} * A'  (m x m)
%   3. Solve Sc * lambda = A * H^{-1} * g - b  for lambda
%   4. Recover x = H^{-1} * (g - A' * lambda)
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

% --- Step 1: Sparse LU factorisation of H ---
% Factorisation satisfies: L * U = R(:,p) * H(:,q)
% Solution to H * x = rhs recovered as: x(q) = U \ (L \ (R(:,p) \ rhs))
[L, U, p, q, R] = lu(H, 'vector');

% --- Step 2: Form Schur complement Sc = A * H^{-1} * A' ---
HinvAt      = sparse(n, m);
HinvAt(q,:) = U \ (L \ (R(:,p) \ A'));
Sc          = A * HinvAt;            % m x m Schur complement

% --- Step 3: Solve for multipliers lambda ---
% A * H^{-1} * A' * lambda = A * H^{-1} * g - b
Hinvg    = zeros(n, 1);
Hinvg(q) = U \ (L \ (R(:,p) \ g));
lambda   = Sc \ (A * Hinvg - b);

% --- Step 4: Recover primal solution x ---
% x = H^{-1} * (g - A' * lambda)
x    = zeros(n, 1);
x(q) = U \ (L \ (R(:,p) \ (g - A' * lambda)));

end




