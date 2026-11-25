function WriteInverseRestartFileDuringInversion(p,Fx,J,dJdp,JGHouts,RunInfo,UserVar,CtrlVar,MUA,BCs,l,InvStartValues,Priors,Meas,BCsAdjoint,Ftemplate)
%%
% Writes an inverse restart file during an inversion, with the current iterate as the final values, so that the inversion can be
% continued from it (eg if the run is stopped). Called by the optimiser at intervals of CtrlVar.Inverse.RestartFileWriteInterval
% iterations (currently -UaOptimisation-HessianVectorBased-, see InvertForModelParameters.m). The variables of the inversion are not
% changed: the state is copied into a local F. (23 November 2026)
%
%   p        : current values of the parameters (the optimisation variables)
%   Fx       : the state returned by the cost function at p
%   J, dJdp  : cost function and gradient at p
%   JGHouts  : additional outputs of the cost function at p (I, R, ...)
%   RunInfo  : current RunInfo, including the iteration history and the trust-region state (for continuing in a restart run)
%   Ftemplate: the state at the start of the inversion, used if Fx is empty (optimisers that do not keep the state at p)
%
if isempty(Fx)
    Fx=Ftemplate;
end
F=p2F(CtrlVar,MUA,p,Fx,Meas,Priors);
InvFinalValues=Vars2InvValues(CtrlVar,F,InvStartValues,J,dJdp,JGHouts,RunInfo,[]);
WriteAdjointRestartFile(UserVar,CtrlVar,MUA,BCs,F,F.GF,l,RunInfo,InvStartValues,Priors,Meas,BCsAdjoint,InvFinalValues);
end
