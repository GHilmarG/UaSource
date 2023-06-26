function  ResidualsCriteria=uvResidualsCriteria(CtrlVar,rForce,rWork,iteration,gamma)


if gamma > max(CtrlVar.uvExitBackTrackingStepLength,CtrlVar.BacktrackingGammaMin)

    % This is the desired tolerance

    ResidualsCriteria=(rWork<CtrlVar.uvDesiredWorkAndForceTolerances(1)  && rForce<CtrlVar.uvDesiredWorkAndForceTolerances(2))...
        && (rWork<CtrlVar.uvDesiredWorkOrForceTolerances(1)  || rForce<CtrlVar.uvDesiredWorkOrForceTolerances(2))...
        && iteration >= CtrlVar.NRitmin;


else

    % This is the acceptable tolerance

    ResidualsCriteria=(rWork<CtrlVar.uvAcceptableWorkAndForceTolerances(1)  && rForce<CtrlVar.uvAcceptableWorkAndForceTolerances(2))...
        && (rWork<CtrlVar.uvAcceptableWorkOrForceTolerances(1)  || rForce<CtrlVar.uvAcceptableWorkOrForceTolerances(2))...
        && iteration >= CtrlVar.NRitmin;

end

end