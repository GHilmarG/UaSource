function [x, lambda] = KKT_null_space_decomposition_scaled(H, g, A, b)
% KKT_NULL_SPACE_SCALED_V2  Null-space solution of the KKT system, with a
%                           fast path for selection-matrix constraints.
%
% Solves
%
%   | H   A' | | x      |   | g |
%   | A   0  | | lambda | = | b |
%
% Fast path
% ---------
% When every row of A has exactly one nonzero - which is the case whenever
% the constraints fix individual degrees of freedom, as thickness-minimum
% constraints do - the null-space basis Z is a selection matrix.  Three
% steps of the general algorithm then collapse:
%
%   * Z need not be formed by QR at all; it is the identity restricted to
%     the unconstrained ("free") rows, so it is just an index list.
%   * A*A' is diagonal, so the factorisation used for the particular
%     solution and for multiplier recovery becomes a diagonal solve.
%   * The reduced Hessian Z'*Hs*Z is simply Hs(free,free), pure indexing
%     rather than two sparse matrix products.
%
% What remains is one factorisation of an (n-m) x (n-m) matrix, which is
% smaller than the n x n system that pre-elimination factorises.
%
% The general path (QR-based) is retained for constraints that couple
% several degrees of freedom.
%
% Scaling
% -------
% Symmetric Jacobi scaling, Dx = diag(1./sqrt(|diag(H)|)), is applied to H,
% and row scaling to A.
%
% This is deliberately NOT the row-2-norm scaling used in
% KKT_null_space_decomposition_scaled.m.  Measured on a uvh system with
% n=330291, m=563:
%
%     condest(H)                 6.0e+12
%     Jacobi   Dx*H*Dx           1.1e+07
%     row-norm Dx*H*Dx           4.8e+13
%
% so row-norm scaling leaves the matrix eight times worse conditioned than
% doing nothing at all, while Jacobi scaling improves it by nearly six
% orders of magnitude.
%
% The reason is that the ill-conditioning does not come from the different
% units of velocity (m/yr) and thickness (m) - the three diagonal blocks
% sit within a factor of 23 of each other.  It comes from within the
% momentum block, whose diagonal spans eleven orders of magnitude because
% h*eta varies that much between thick fast-flowing and thin nearly
% stagnant ice.  Scaling by the diagonal targets exactly that; scaling by
% row norms does not, because a row norm is dominated by the largest
% off-diagonal entry rather than by the local magnitude of h*eta.
%
% Note that with a direct solver this changes accuracy very little, since
% UMFPACK already applies its own row scaling internally - the measured
% error against a refined solution was 4.7e-13 unscaled against 2.7e-13
% scaled.  The scaling matters for the conditioning of the reduced system
% and for any iterative method applied to it.
%
% Inputs / outputs as KKT_null_space_decomposition_scaled.

[m, n] = size(A);
assert(n == size(H,1),  'H and A must have compatible dimensions');
assert(n == size(H,2),  'H must be square');
assert(m == length(b),  'A and b must have compatible dimensions');
assert(n == length(g),  'H and g must have compatible dimensions');
assert(m <= n,          'System must not be over-determined (m <= n)');

% --- Step 1: scaling -------------------------------------------------
% Jacobi scaling on H: target the diagonal, which is where the eleven
% orders of magnitude of spread live.  Fall back to the row 2-norm for any
% row whose diagonal entry is zero or non-finite, and to unity if that
% fails too, so the scaling can never introduce Inf or NaN.
dH = full(abs(diag(H)));
bad = ~isfinite(dH) | dH <= 0;
if any(bad)
    rn = full(sqrt(sum(H(bad,:).^2,2)));
    rn(~isfinite(rn) | rn <= 0) = 1;
    dH(bad) = rn;
end
dx = 1./sqrt(dH);

% Row scaling on A.  For a selection matrix this simply normalises each
% constraint row to unit magnitude.
dc = full(sqrt(sum(A.^2,2)));
dc(~isfinite(dc) | dc <= 0) = 1;
dc = 1./dc;

Dx = spdiags(dx,0,n,n);
Dc = spdiags(dc,0,m,m);

Hs = Dx*H*Dx;
As = Dc*A*Dx;
gs = Dx*g;
bs = Dc*b;

isSelection = all(sum(A~=0,2)==1);

if isSelection

    % ---- fast path --------------------------------------------------
    [rowA,colA] = find(As);
    [~,ord]     = sort(rowA);
    cons        = colA(ord);            % constrained column, per constraint row
    aval        = full(As(sub2ind([m n],(1:m)',cons)));

    free            = true(n,1);
    free(cons)      = false;

    % particular solution: put the constrained values in place, zero elsewhere
    xps        = zeros(n,1);
    xps(cons)  = bs./aval;

    % reduced system: delete constrained rows and columns
    Hr   = Hs(free,free);
    rhs  = gs(free) - Hs(free,:)*xps;

    dHr  = decomposition(Hr,'auto');
    u    = dHr\rhs;

    xs        = xps;
    xs(free)  = xs(free) + u;

    % multipliers: As*As' is diagonal, so this is a scalar divide per row
    r        = gs - Hs*xs;
    lambdas  = r(cons)./aval;

else

    % ---- general path -----------------------------------------------
    [Q,~] = qr(As');
    Z     = Q(:,m+1:end);

    dAAsT = decomposition(As*As','auto');
    xps   = As'*(dAAsT\bs);

    Hr  = Z'*(Hs*Z);
    dHr = decomposition(Hr,'auto');

    u  = dHr\(Z'*(gs - Hs*xps));
    xs = xps + Z*u;

    lambdas = dAAsT\(As*(gs - Hs*xs));

end

% --- unscale ----------------------------------------------------------
x      = Dx*xs;
lambda = Dc*lambdas;

end
