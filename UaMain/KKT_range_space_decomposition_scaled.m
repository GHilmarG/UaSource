





function [x, lambda] = KKT_range_space_decomposition_scaled(H, g, A, b)
% KKT_RANGE_SPACE_SCALED  Solves the KKT system using the range space method
%                         with diagonal scaling and the MATLAB decomposition
%                         object.
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
%   2. Factorise Hs using the MATLAB decomposition object
%   3. Form scaled Schur complement Sc = As * Hs^{-1} * As'  (m x m)
%   4. Solve Sc * lambdas = As * Hs^{-1} * gs - bs  for lambdas
%   5. Recover xs = Hs^{-1} * (gs - As' * lambdas)
%   6. Unscale to recover x and lambda
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


% --- Step 2: Factorise scaled Hs ---
dH = decomposition(Hs);

% --- Step 3: Form scaled Schur complement Sc = As * Hs^{-1} * As' ---
HinvAt = dH \ (As');             % n x m
Sc     = As * HinvAt;                % m x m



% --- Step 4: Solve for scaled multipliers ---
% As * Hs^{-1} * As' * lambdas = As * Hs^{-1} * gs - bs
Hinvgs  = dH \ gs;
lambdas = Sc \ (As * Hinvgs - bs);

% --- Step 5: Recover scaled primal solution ---
% xs = Hs^{-1} * (gs - As' * lambdas)
xs = dH \ (gs - As' * lambdas);

% --- Step 6: Unscale ---
x      = Dx * xs;
lambda = Dc * lambdas;

% --- Diagnostic residuals (comment out in production) ---
res1 = norm(H*x + A'*lambda - g)/norm(x);
res2 = norm(A*x - b)/norm(x);
fprintf('Stationarity  ||H*x + A''*lambda - g||/||x||= %.2e\n', res1);
fprintf('Feasibility   ||A*x - b||/||x||             = %.2e\n', res2);


end