function  ResidualsCriteria=uvResidualsCriteria(CtrlVar,rForce,rWork,iteration,gamma) %#ok<INUSL>

% Convergence test of the uv solve (SSTREAM2dNR2). The desired tolerance is used after a normal step, the acceptable tolerance if the
% step was very short (gamma below CtrlVar.uvExitBackTrackingStepLength), or if the calling routine passes gamma=0 (failed backtracking,
% stagnation), as in the uvh solver (uvhResidualsCriteria.m).
%
% (1 Nov 2026) The criterion depends on rForce only (CtrlVar.uvDesiredForceTolerance, CtrlVar.uvAcceptableForceTolerance). Previously, rWork
% also entered (CtrlVar.uv*WorkAndForceTolerances, CtrlVar.uv*WorkOrForceTolerances, with all work tolerances inf, so that only a NaN rWork,
% returned by the line search after a step that is not a Newton step, could change the decision: NaN<inf is false, so that the iteration
% could not be judged converged after such a step). rWork is kept in the argument list so that the calls need not change.

if gamma > max(CtrlVar.uvExitBackTrackingStepLength,CtrlVar.BacktrackingGammaMin)
    ResidualsCriteria = rForce < CtrlVar.uvDesiredForceTolerance && iteration >= CtrlVar.NRitmin;      % desired tolerance
else
    ResidualsCriteria = rForce < CtrlVar.uvAcceptableForceTolerance && iteration >= CtrlVar.NRitmin;   % acceptable tolerance
end

end
