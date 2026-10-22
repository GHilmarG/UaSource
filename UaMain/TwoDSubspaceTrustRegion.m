



function [p,y,info] = TwoDSubspaceTrustRegion(g,H,v1,v2,Delta)




%%
%   Exact solution of the trust-region subproblem restricted
%   to the 2D subspace spanned by v1 and v2.
%
%   min_p   m(p) = g'*p + (1/2) p'*H*p     s.t.  ||p|| <= Delta
%   p in span{v1,v2}
%
%   [p,y,info] = TwoDSubspaceTR(g,H,v1,v2,Delta)
%
%   INPUTS
%     g     : n-by-1 gradient
%     H     : n-by-n Hessian (numeric matrix) OR a function handle
%             Hfun(x) returning H*x  (Hessian-vector product)
%     v1,v2 : n-by-1 vectors spanning the subspace (need not be
%             orthogonal or normalised, e.g. v1=-g, v2=Newton step)
%     Delta : trust-region radius (||p|| <= Delta, Euclidean norm)
%
%   OUTPUTS
%     p     : n-by-1 solution in the full space, p = y1*v1 + y2*v2
%     y     : 2-by-1 coefficients, p = V*y, V=[v1 v2]
%     info  : struct with diagnostic fields (case, lambda, mu, etc.)
%
%   Method: Byrd, Schnabel & Shultz (1988) two-dimensional subspace
%   minimisation. The (generally non-orthogonal) basis is whitened using
%   the Gram matrix so that the norm constraint becomes a standard
%   Euclidean ball, the reduced 2x2 Hessian is diagonalised, and the
%   boundary solution is obtained as the root of a quartic secular
%   equation. Candidates are selected by direct evaluation of the
%   reduced quadratic model, which is more robust than trying to reason
%   about "the" correct root analytically.
%%


tol = 1e-10;

% ---- basic setup -----------------------------------------------
v1 = v1(:); v2 = v2(:);
V  = [v1, v2];

if isa(H,'function_handle')
    Hv1 = H(v1);
    Hv2 = H(v2);
else
    Hv1 = H*v1;
    Hv2 = H*v2;
end

Ghat = V'*g(:);                 % 2x1 reduced gradient
Hhat = [v1'*Hv1, v1'*Hv2; ...
        v2'*Hv1, v2'*Hv2];      % 2x2 reduced Hessian
Hhat = (Hhat+Hhat')/2;          % enforce symmetry (roundoff)

M = V'*V;                       % 2x2 Gram matrix (metric for ||p||)
M = (M+M')/2;

info = struct('case','','lambda',0,'mu',[],'a',[]);

% Reduced cost function value (up to the constant f), for candidate
% comparison. y is 2x1.
mval = @(y) Ghat'*y + 0.5*y'*Hhat*y;

% ---- 1. Try the interior (unconstrained) minimiser -------------
% Only valid if Hhat is SPD and the resulting point lies inside Delta.
[~,pdflag] = chol(Hhat,'lower');
if pdflag==0
    yU = -(Hhat\Ghat);
    if yU'*M*yU <= Delta^2*(1+1e-10)
        y = yU;
        p = V*y;
        info.case = 'interior';
        return
    end
end

% ---- 2. Boundary solution: whiten the constraint ----------------
% M = L*L', z = L'*y  =>  ||p||^2 = y'*M*y = z'*z
L = chol(M,'lower');            % M is SPD as long as v1,v2 independent

Hz = L\(Hhat/L');               % = L^{-1} Hhat L^{-T}
Hz = (Hz+Hz')/2;
gz = L\Ghat;

% ---- 3. Diagonalise the whitened reduced Hessian -----------------
[Q,Lam] = eig(Hz);
mu = diag(Lam);
[mu,idx] = sort(mu,'ascend');   % mu(1) <= mu(2)
Q = Q(:,idx);
a = Q'*gz;                      % a(1), a(2)

mu1 = mu(1); mu2 = mu(2);
a1  = a(1);  a2  = a(2);

info.mu = mu; info.a = a;

candidates = {};  % each entry: y (2x1) candidate in ORIGINAL reduced coords

% ---- 4. Hard case check ------------------------------------------
% Hard case: most negative-curvature direction has (numerically) no
% gradient component, and mu1<0, so lambda=-mu1 solves secular eqn
% for any component along that eigendirection.
gscale = norm(a) + eps;
if mu1 < 0 && abs(a1) <= tol*max(1,gscale)
    lambda = -mu1;
    info.case = 'hard case';
    info.lambda = lambda;

    denom2 = mu2+lambda;
    if abs(denom2) > tol
        w2 = -a2/denom2;
    else
        w2 = 0; % fully degenerate (mu1=mu2), any direction works
    end
    w1sq = max(Delta^2 - w2^2,0);
    w1 = sqrt(w1sq);

    for s = [1,-1]
        w = [s*w1; w2];
        y_cand = L'\(Q*w);
        candidates{end+1} = y_cand; %#ok<AGROW>
    end
else
    % ---- 5. Standard boundary case: solve the quartic -------------
    info.case = 'boundary (quartic)';

    s = mu1+mu2;
    q = mu1*mu2;

    c4 = Delta^2;
    c3 = 2*Delta^2*s;
    c2 = Delta^2*(s^2+2*q) - (a1^2+a2^2);
    c1 = 2*Delta^2*s*q - 2*(a1^2*mu2 + a2^2*mu1);
    c0 = Delta^2*q^2 - (a1^2*mu2^2 + a2^2*mu1^2);

    r = roots([c4 c3 c2 c1 c0]);

    lamMin = max(0,-mu1);  % feasibility requirement: Hz+lambda*I >= 0

    for k = 1:numel(r)
        lambda = r(k);
        if abs(imag(lambda)) > 1e-8*max(1,abs(real(lambda)))
            continue % discard complex roots
        end
        lambda = real(lambda);
        if lambda < lamMin - 1e-8
            continue % infeasible root
        end

        d1 = mu1+lambda; d2 = mu2+lambda;
        if abs(d1) < tol || abs(d2) < tol
            continue % (near) division by zero, skip -- hard case path handles mu1+lambda=0
        end
        w = [-a1/d1; -a2/d2];

        y_cand = L'\(Q*w);
        candidates{end+1} = y_cand; 
    end

    if isempty(candidates)
        % Fallback: numerical trouble with the quartic. Use the
        % steepest-descent-like boundary point in the whitened space
        % along -gz direction, scaled to the boundary.
        if norm(gz) > 0
            zfallback = -Delta*gz/norm(gz);
        else
            zfallback = [Delta;0];
        end
        y_cand = L'\zfallback;
        candidates{end+1} = y_cand; 
        info.case = [info.case ' (fallback)'];
    end
end

% ---- 6. Pick the best candidate by direct evaluation -------------
bestVal = inf; bestY = candidates{1};
for k = 1:numel(candidates)
    yk = candidates{k};
    % re-normalise onto the exact boundary to clean up roundoff
    nrm = sqrt(yk'*M*yk);
    if nrm > 0
        yk = yk*(Delta/nrm);
    end
    val = mval(yk);
    if val < bestVal
        bestVal = val;
        bestY = yk;
    end
end

y = bestY;
p = V*y;
end