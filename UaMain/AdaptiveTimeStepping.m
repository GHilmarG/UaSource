function [RunInfo,dtOut,dtRatio]=AdaptiveTimeStepping(UserVar,RunInfo,CtrlVar,MUA,F)

%% [RunInfo,dtOut,dtRatio]=AdaptiveTimeStepping(UserVar,RunInfo,CtrlVar,MUA,F)
%
% Automated time stepping. The time step on input is F.dt, and dtOut is the new time step. It is up to the calling program to set F.dt=dtOut.
%
% (10 Oct 2026) This is a driver that selects between the available approaches, depending on CtrlVar.AdaptiveTimeSteppingMethod:
%
%   "iteration-based" : the time step is based on the number of non-linear iterations (AdaptiveTimeSteppingIterationBased.m). This was
%                       the approach used in AdaptiveTimeStepping.m previously, and is the default.
%   "error-estimate"  : the time step is based on an estimate of the time-discretisation error of the ice thickness, together with the
%                       number of non-linear iterations (AdaptiveTimeSteppingWithErrorEstimate.m).
%
%%

narginchk(5,5)

Method="iteration-based";
if isfield(CtrlVar,"AdaptiveTimeSteppingMethod")
    Method=string(CtrlVar.AdaptiveTimeSteppingMethod);
end

switch Method

    case "iteration-based"
        [RunInfo,dtOut,dtRatio]=AdaptiveTimeSteppingIterationBased(UserVar,RunInfo,CtrlVar,MUA,F);

    case "error-estimate"
        [RunInfo,dtOut,dtRatio]=AdaptiveTimeSteppingWithErrorEstimate(UserVar,RunInfo,CtrlVar,MUA,F);

    otherwise
        error("AdaptiveTimeStepping:UnknownMethod","Unknown value of CtrlVar.AdaptiveTimeSteppingMethod: %s. Use ""iteration-based"" or ""error-estimate"".",Method)

end

end
