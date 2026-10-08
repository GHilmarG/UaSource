function ResidualsCriteria=uvhResidualsCriteria(CtrlVar,rForce,rWork,iteration,Acceptable)
%%
% ResidualsCriteria=uvhResidualsCriteria(CtrlVar,rForce,rWork,iteration,Acceptable)
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
%
% (8 Oct 2026) Moved here from SSTREAM_TransientImplicit.
%
%%

narginchk(5,5)

if isnan(rWork)
    rWork=0;
end

if ~Acceptable

    ResidualsCriteria=(rWork<CtrlVar.uvhDesiredWorkAndForceTolerances(1)  && rForce<CtrlVar.uvhDesiredWorkAndForceTolerances(2))...
        && (rWork<CtrlVar.uvhDesiredWorkOrForceTolerances(1)  || rForce<CtrlVar.uvhDesiredWorkOrForceTolerances(2))...
        && iteration >= CtrlVar.NRitmin;

else

    ResidualsCriteria=(rWork<CtrlVar.uvhAcceptableWorkAndForceTolerances(1)  && rForce<CtrlVar.uvhAcceptableWorkAndForceTolerances(2))...
        && (rWork<CtrlVar.uvhAcceptableWorkOrForceTolerances(1)  || rForce<CtrlVar.uvhAcceptableWorkOrForceTolerances(2))...
        && iteration >= CtrlVar.NRitmin;

end

end
