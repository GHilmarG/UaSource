function [RunInfo,dtAccuracy]=ATSTimeDiscretisationErrorLimit(CtrlVar,RunInfo)
%%
% [RunInfo,dtAccuracy]=ATSTimeDiscretisationErrorLimit(CtrlVar,RunInfo)
%
% Accuracy-based upper limit on the next time step, based on the estimated local time-discretisation error of the ice thickness in the
% last time step (see TimeDiscretisationErrorEstimate.m). Called from AdaptiveTimeStepping.m, where the new time step is the minimum of
% this limit and the time step proposed by the automated time stepping based on the number of non-linear iterations. (10 Oct 2026)
%
% The error measure is the error per unit time (EPUS), scaled by the tolerances:
%
%   q = area-weighted rms over the nodes where the estimate is defined of (e/dt)./(atol+rtol*|h|)
%
% where e is the estimated local error of h in the last time step, dt the length of that time step, and atol (m/yr) and rtol (1/yr) are given
% in CtrlVar.ATSTimeDiscretisationError. q<=1 means that the tolerance is met. For the second-order time stepping (theta=0.5) the local error
% per unit time is proportional to dt^2, and the time step for which q=1 is dt*q^(-1/2). The limit is
%
%   dtAccuracy = Safety * dt * q^(-1/2)
%
% bounded by MaxIncrease*dt from above and dt/MaxDecrease from below. Integrated over a run of length T, the accumulated local error is then
% approximately bounded by T*atol (for rtol=0).
%
% No limit (dtAccuracy=Inf) is returned, ie the time step is then determined by the other criteria alone, if in the last time step:
%   - no estimate was calculated (eg in the first two time steps, or for theta~=0.5)
%   - the area where the estimate is defined was smaller than the fraction MinValidArea of the total area
%   - the time was less than StartTime (eg to exclude an initial adjustment period)
%
% The values of q and dtAccuracy are stored in RunInfo.Forward.ATSTimeDiscretisationErrorRatio and RunInfo.Forward.ATSTimeDiscretisationErrorDtLimit,
% indexed by the run-step number of the time step for which the limit is calculated.
%
%%

dtAccuracy=Inf;

Opt=CtrlVar.ATSTimeDiscretisationError;
Def=struct("atol",0.01,"rtol",0,"Safety",0.85,"MaxIncrease",2,"MaxDecrease",5,"StartTime",-inf,"MinValidArea",0.1,"InfoLevel",1);
fn=fieldnames(Def);
for I=1:numel(fn)
    if ~isfield(Opt,fn{I}) ; Opt.(fn{I})=Def.(fn{I}) ; end
end

kNew=CtrlVar.CurrentRunStepNumber;     % the time step for which the time step is now being determined
kLast=kNew-1;                          % the last time step
q=NaN;
Reason="";

Fw=RunInfo.Forward;
Needed=["hTimeDiscretisationErrorFlag","hTimeDiscretisationErrorRateScaled","hTimeDiscretisationErrorDt","hTimeDiscretisationErrorValidArea","hTimeDiscretisationErrorTime"];
if kLast<1 || ~all(isfield(Fw,Needed)) || any(arrayfun(@(n) numel(Fw.(n))<kLast,Needed))
    Reason="no estimate available for the last time step";
else
    Flag=Fw.hTimeDiscretisationErrorFlag(kLast);
    q=Fw.hTimeDiscretisationErrorRateScaled(kLast);
    dtLast=Fw.hTimeDiscretisationErrorDt(kLast);
    ValidArea=Fw.hTimeDiscretisationErrorValidArea(kLast);
    tLast=Fw.hTimeDiscretisationErrorTime(kLast);

    if Flag~=0 || ~isfinite(q) || ~(dtLast>0)
        Reason="no estimate for the last time step";
    elseif ValidArea<Opt.MinValidArea
        Reason=sprintf("estimate defined over only %.1f%% of the area",100*ValidArea);
    elseif tLast<Opt.StartTime
        Reason=sprintf("time is less than StartTime=%g",Opt.StartTime);
    else
        if q>0
            Factor=Opt.Safety*q^(-1/2);
        else
            Factor=Opt.MaxIncrease;
        end
        Factor=min(max(Factor,1/Opt.MaxDecrease),Opt.MaxIncrease);
        dtAccuracy=Factor*dtLast;
    end
end

RunInfo=SetLimitSeries(RunInfo,kNew,q,dtAccuracy);

if Opt.InfoLevel>=1
    if isfinite(dtAccuracy)
        fprintf(' Accuracy-based time step: scaled time-discretisation error per unit time q=%-9.3g (q<=1: tolerance met), time step allowed by accuracy %-g \n',q,dtAccuracy)
    else
        fprintf(' Accuracy-based time step: no limit (%s). \n',Reason)
    end
end

end


function RunInfo=SetLimitSeries(RunInfo,k,q,dtAccuracy)
% Stores q and dtAccuracy for run step k in RunInfo.Forward, creating or extending the series with NaN if needed.
Names=["ATSTimeDiscretisationErrorRatio","ATSTimeDiscretisationErrorDtLimit"];
Values=[q,dtAccuracy];
for I=1:numel(Names)
    if ~isfield(RunInfo.Forward,Names(I))
        RunInfo.Forward.(Names(I))=NaN(max(k,numel(RunInfo.Forward.time)),1);
    elseif numel(RunInfo.Forward.(Names(I)))<k
        RunInfo.Forward.(Names(I))=[RunInfo.Forward.(Names(I))(:);NaN(k-numel(RunInfo.Forward.(Names(I))),1)];
    end
    RunInfo.Forward.(Names(I))(k)=Values(I);
end
end
