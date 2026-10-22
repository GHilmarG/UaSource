

function [p0_new, Delta_new, accepted] = TrustRegionUpdate(p0, dp, J0, Jnew, pred, Delta, DeltaMax,stepnorm)


%%
% Standard trust-region ratio test and radius update based on Nocedal & Wright and various others
%   p0     : current point
%   dp      : proposed step (from TwoDSubspaceTR)
%   J0     : J(p0)            (already have this)
%   Jnew   : J(p0+p)          (one extra function evaluation)
%   pred   : predicted reduction, -(g'*p + 0.5*p'*H*p), should be >= 0
%   Delta  : current trust-region radius
%   DeltaMax : cap on radius growth
%%


eta1 = 0.25;   % shrink threshold
eta2 = 0.75;   % grow threshold
eta0 = 1e-4;   % minimum acceptance threshold (accept step if rho > eta0)
gamma1 = 0.25; % shrink factor
gamma2 = 2.0;  % grow factor

ared = J0 - Jnew;

if pred <= 0
    % Should not happen for a correctly computed subspace step;
    % treat conservatively as a rejected, model-breakdown step.
    rho = -inf;
else
    rho = ared / pred;
end

%stepnorm = norm(dp);


% Update Radius
if rho < eta1
    Delta_new = gamma1 * Delta;
elseif rho > eta2 && stepnorm >= (1-1e-8)*Delta   % step was at the boundary
    Delta_new = min(gamma2*Delta, DeltaMax);
else
    Delta_new = Delta;                           
end

% accept/reject the step itself: 
if rho > eta0
    p0_new = p0 + dp;
    accepted = true;
else
    p0_new = p0;                                 
    accepted = false;
end
end