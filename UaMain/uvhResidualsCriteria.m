function ResidualsCriteria=uvhResidualsCriteria(CtrlVar,rForce,rWork,iteration,Acceptable,rForceFloor)
%%
% ResidualsCriteria=uvhResidualsCriteria(CtrlVar,rForce,rWork,iteration,Acceptable,rForceFloor)
%
% Convergence test for the uvh system. Used both by the implicit uvh solver (SSTREAM_TransientImplicit) and by the semi-implicit
% uv-h solver (uvhSemiImplicit), so that both solvers use exactly the same convergence criterion.
%
%   rForce, rWork : as returned by CalcCostFunctionNRuvh
%   iteration     : (outer) iteration number, compared with CtrlVar.NRitmin
%   Acceptable    : if true, the 'acceptable' tolerances are used, otherwise the 'desired' tolerances:
%
%       CtrlVar.uvhDesiredWorkAndForceTolerances,    CtrlVar.uvhDesiredWorkOrForceTolerances
%       CtrlVar.uvhAcceptableWorkAndForceTolerances, CtrlVar.uvhAcceptableWorkOrForceTolerances
%
%   In SSTREAM_TransientImplicit the acceptable tolerances are used if the last backtracking step was short.
%   If rWork is nan (rWork has no meaning for some minimisation approaches) it is not used, ie set to zero.
%   rForceFloor   : (optional, default 0) estimate of the round-off floor of rForce (see uvhResidualRoundOffFloor). The force
%                   tolerances are replaced by max(tolerance,rForceFloor), so that convergence can be reached for small time steps.
%
% (8 Oct 2026) Moved here from SSTREAM_TransientImplicit.
%
%%

narginchk(5,6)

if nargin<6 || isempty(rForceFloor)
    rForceFloor=0;
end

if isnan(rWork)
    rWork=0;
end

if ~Acceptable

    ResidualsCriteria=(rWork<CtrlVar.uvhDesiredWorkAndForceTolerances(1)  && rForce<max(CtrlVar.uvhDesiredWorkAndForceTolerances(2),rForceFloor))...
        && (rWork<CtrlVar.uvhDesiredWorkOrForceTolerances(1)  || rForce<max(CtrlVar.uvhDesiredWorkOrForceTolerances(2),rForceFloor))...
        && iteration >= CtrlVar.NRitmin;

else

    ResidualsCriteria=(rWork<CtrlVar.uvhAcceptableWorkAndForceTolerances(1)  && rForce<max(CtrlVar.uvhAcceptableWorkAndForceTolerances(2),rForceFloor))...
        && (rWork<CtrlVar.uvhAcceptableWorkOrForceTolerances(1)  || rForce<max(CtrlVar.uvhAcceptableWorkOrForceTolerances(2),rForceFloor))...
        && iteration >= CtrlVar.NRitmin;

end

end
