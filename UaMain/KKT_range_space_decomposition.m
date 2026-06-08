
function [x, lambda] = KKT_range_space_decomposition(H, g, A, b)
% KKT_RANGE_SPACE_DECOMP  Solves the KKT system using the range space method
%                         with the MATLAB decomposition object.
%
% Solves the KKT system:
%
%   | H   A' | | x      |   | g |
%   | A   0  | | lambda | = | b |
%
% Method:
%   1. Factorise H using the MATLAB decomposition object
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

% --- Step 1: Factorise H ---

dH = decomposition(H);

% --- Step 2: Form Schur complement Sc = A * H^{-1} * A' ---
% Solve H * X = A' for all m right hand sides simultaneously
HinvAt = dH \ (A');                  % n x m:  H^{-1} * A'
Sc     = A * HinvAt;                 % m x m:  A * H^{-1} * A'

% --- Step 3: Solve for multipliers lambda ---
% Substituting x = H^{-1}(g - A'*lambda) into A*x = b gives:
% A * H^{-1} * A' * lambda = A * H^{-1} * g - b
Hinvg  = dH \ g;                     % n x 1:  H^{-1} * g
lambda = Sc \ (A * Hinvg - b);       % m x 1

% --- Step 4: Recover primal solution ---
% x = H^{-1} * (g - A' * lambda)
x = dH \ (g - A' * lambda);          % n x 1

% --- Diagnostic residuals (comment out in production) ---
res1 = norm(H*x + A'*lambda - g);
res2 = norm(A*x - b);
fprintf('Stationarity  ||H*x + A''*lambda - g|| = %.2e\n', res1);
fprintf('Feasibility   ||A*x - b||              = %.2e\n', res2);


end