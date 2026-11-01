function ResidualsCriteria=uvhResidualsCriteria(CtrlVar,rForce,rWork,iteration,Acceptable,rForceFloor) %#ok<INUSL>
%%
% ResidualsCriteria=uvhResidualsCriteria(CtrlVar,rForce,rWork,iteration,Acceptable,rForceFloor)
%
% Convergence test for the uvh system. Used both by the implicit uvh solver (SSTREAM_TransientImplicit) and by the semi-implicit
% uv-h solver (uvhSemiImplicit), so that both solvers use exactly the same convergence criterion.
%
%   rForce        : as returned by CalcCostFunctionNRuvh
%   rWork         : not used (see below); kept in the argument list so that the calls need not change
%   iteration     : (outer) iteration number, compared with CtrlVar.NRitmin
%   Acceptable    : if true, the 'acceptable' tolerance is used, otherwise the 'desired' tolerance:
%
%       CtrlVar.uvhDesiredForceTolerance     (default 1e-15)
%       CtrlVar.uvhAcceptableForceTolerance  (default 1e-6)
%
%   In SSTREAM_TransientImplicit the acceptable tolerance is used if the iteration has stagnated or the last backtracking step was
%   short.
%   rForceFloor   : (optional, default 0) estimate of the round-off floor of rForce (see uvhResidualRoundOffFloor). The tolerances are
%                   replaced by max(tolerance,rForceFloor), so that convergence can be reached for small time steps.
%
% (8 Oct 2026) Moved here from SSTREAM_TransientImplicit.
%
% (1 Nov 2026) The criterion depends on rForce only. Previously, rWork also entered (CtrlVar.uvhAcceptableWorkOrForceTolerances=[1 1e-8]
% required rWork<1 or rForce<1e-8 for acceptance, in addition to rForce<1e-6). rWork contains the Newton step, which becomes very large
% at nearly singular modes of the Jacobian (eg weakly determined velocities of floating, nearly ice-free areas), and is returned as NaN
% (and then set to zero) by the line search whenever the accepted step is not a Newton step. The acceptance of a stalled solve therefore
% depended on which line-search branch produced the last step (Thule: a solve stalled at rForce=3e-7 after a Newton step was rejected,
% one stalled at rForce=5e-7 after a fallback step was accepted). The previous fields CtrlVar.uvh*WorkAndForceTolerances and
% CtrlVar.uvh*WorkOrForceTolerances are no longer used (see CtrlVarValidityCheck.m).
%
%%

narginchk(5,6)

if nargin<6 || isempty(rForceFloor)
    rForceFloor=0;
end

if ~Acceptable
    ResidualsCriteria = rForce < max(CtrlVar.uvhDesiredForceTolerance,rForceFloor) && iteration >= CtrlVar.NRitmin ;
else
    ResidualsCriteria = rForce < max(CtrlVar.uvhAcceptableForceTolerance,rForceFloor) && iteration >= CtrlVar.NRitmin ;
end

end
