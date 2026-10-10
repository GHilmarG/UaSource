function [RunInfo,F1]=TimeDiscretisationErrorEstimate(CtrlVar,RunInfo,MUA,F0,F1,Fm1,BCs1,hPosNodeStartOfStep)
%%
% [RunInfo,F1]=TimeDiscretisationErrorEstimate(CtrlVar,RunInfo,MUA,F0,F1,Fm1,BCs1,hPosNodeStartOfStep)
%
% Estimate of the local time-discretisation error of the ice thickness in the last time step. (10 Oct 2026)
%
% On return F1.hTimeDiscretisationErrorEstimate is a nodal field with the estimated local error (in units of h, ie m) of the time step
% from t(n) to t(n+1)=t(n)+dt, ie the estimated difference between the thickness obtained by an exact time integration over this step
% (starting from the computed thickness at t(n)) and the computed thickness at t(n+1). Norms of this field are stored as time series in
% RunInfo.Forward (see below).
%
% Method: Predictor-corrector estimate (Milne's device). The predictor is the second-order explicit estimate of h at t(n+1) obtained from
% backward differences over the two previous time steps (see ExplicitEstimationUsingBackwardDifferences.m), and the corrector is the
% computed (implicit) solution using the trapezoidal rule (theta=0.5). For a solution h(t) that is smooth in time:
%
%   h(t(n+1)) - hC  ~  CC h''',  CC=-dt^3/12                                   (local error of the trapezoidal rule)
%   h(t(n+1)) - hP  ~  CP h''' + W h''',  CP=dt (dt+dtm1) (dt+dtm1+dtm2)/6     (error of the quadratic extrapolation)
%
% where dtm1 and dtm2 are the two previous time steps, and W h''' is the contribution of the local errors of these two previous time
% steps to the extrapolated value (the predictor extrapolates the computed, and not the exact, solution):
%
%   W = w1*c1 + w2*(c1+c2),  c1=-dtm1^3/12, c2=-dtm2^3/12
%
% with w1 and w2 the weights of h(t(n-1)) and h(t(n-2)) in the quadratic extrapolation to t(n+1). Eliminating h''' gives
%
%   e = CC (hC-hP) / (CP-CC-W)        (for a constant time step: e = -(hC-hP)/12)
%
% This was verified on a nonlinear test problem with an exact solution, for both constant and variable time steps: the estimate converges
% to the true local error, and for coarse time steps it is somewhat too large (conservative).
%
% Interpretation: This is an estimate of the local error, ie the error committed in this time step, in metres. The global error at a
% given time is the accumulation of the local errors (some of which may decay), and is typically bounded by the sum of the local errors.
% The estimate assumes that h is smooth in time over the three time steps involved, and is only an indicator where this is not the case,
% for example where the active set changes, or at jumps in the forcing. It does not include spatial-discretisation errors.
%
% The estimate is set to NaN at:
%   - nodes with thickness boundary conditions
%   - nodes in the active set of the thickness constraints, at the start or the end of the time step
%   - nodes whose constraint state changed (activated or released) within the last three time steps, and their neighbours (ie all nodes
%     of elements containing such nodes). (10 Oct 2026) Previously the neighbours of all active nodes were excluded, which excluded
%     for example fast-flowing ice next to ice-free (constrained) areas permanently. Next to nodes that remain constrained the evolution
%     is smooth, and the estimate is valid. Three time steps are used as the predictor uses h at t(n-2), t(n-1) and t(n). The previous
%     rule can be selected with CtrlVar.TimeDiscretisationErrorEstimate.MaskBand="active". As long as the history of changes does not
%     cover three time steps (start of a run, after remeshing), the previous rule is used.
%   - nodes where the thickness, at the start or end of the time step, is at or below ThickMin
%   - nodes within the range of the thickness penalty (if used), ie h < ThickMin+2*delta
%   - nodes outside of the ice front (if the level-set method is used)
%   - nodes where no second-order predictor is available (eg after remeshing).
% It is not calculated at all (empty field) for theta~=0.5 or if the time steps over which the rates were calculated are not available
% (the first time steps of a run, or after a restart from an older restart file).
%
% Time series stored in RunInfo.Forward, indexed by the run-step number:
%
%   hTimeDiscretisationErrorTime      : time at the end of the time step
%   (The norms below are calculated as sqrt(x'*M*x/A), with M the consistent mass matrix, x set to zero at masked nodes, and A the area
%   where the estimate is defined.)
%   hTimeDiscretisationErrorRMS       : rms of the estimate (m)
%   hTimeDiscretisationErrorMax       : maximum absolute value of the estimate (m)
%   hTimeDiscretisationErrorMaxX, ..MaxY : location of the maximum
%   hTimeDiscretisationErrorVolume    : area integral of the estimate (m^3), ie an estimate of the error in ice volume
%   hTimeDiscretisationErrorScaled    : area-weighted rms of e./(atol+rtol*|h|), with atol and rtol in CtrlVar.TimeDiscretisationErrorEstimate
%   hTimeDiscretisationErrorRMSperUnitTime : hTimeDiscretisationErrorRMS/dt
%   hTimeDiscretisationErrorValidArea : fraction of the area where the estimate is defined
%   hTimeDiscretisationErrorFlag      : 0: estimate calculated, 1: theta~=0.5, 2: no second-order predictor, 3: inconsistent array sizes
%   hTimeDiscretisationErrorRateScaled: area-weighted rms of (e/dt)./(atol+rtol*|h|), with atol and rtol (per unit time) in
%                                       CtrlVar.ATSTimeDiscretisationError, used for the accuracy-based time-step limit (ATSTimeDiscretisationErrorLimit.m)
%   hTimeDiscretisationErrorDt        : the time step (the actual length of the time step for which the estimate was made)
%
% Switched on by CtrlVar.TimeDiscretisationErrorEstimate.Use=true.
%
% Accumulated estimate (if CtrlVar.TimeDiscretisationErrorEstimate.Accumulate=true): F1.hTimeDiscretisationErrorAccumulated is a structure with
%
%   Signed    : sum over time steps of the estimate at each node (NaN, ie masked, values are not included)
%   Abs       : sum over time steps of the absolute value of the estimate at each node
%   ValidTime : sum of the time steps for which the estimate was defined at each node
%   TotalTime : total time over which the estimate has been accumulated (also including time steps for which no estimate was calculated)
%
% ValidTime./TotalTime is the fraction of the time covered at each node. Neither Signed nor Abs is the global error, as errors are
% transported and may partly decay, but Signed approximates the global error where the transport is small, and Abs is an upper bound on
% the locally accumulated error. The accumulated estimate is part of F, ie it is saved in restart files and continued after a restart
% (unless CtrlVar.TimeDiscretisationErrorEstimate.ResetAccumulatedAtRestart=true), and it is mapped onto a new mesh after remeshing
% (see MapFbetweenMeshes.m).
%
%%

narginchk(8,8)

F1.hTimeDiscretisationErrorEstimate=[];

if ~(isfield(CtrlVar,"TimeDiscretisationErrorEstimate") && isfield(CtrlVar.TimeDiscretisationErrorEstimate,"Use") && CtrlVar.TimeDiscretisationErrorEstimate.Use)
    return
end

Opt=CtrlVar.TimeDiscretisationErrorEstimate;
if ~isfield(Opt,"atol") ; Opt.atol=1 ; end
if ~isfield(Opt,"rtol") ; Opt.rtol=1e-3 ; end
if ~isfield(Opt,"InfoLevel") ; Opt.InfoLevel=1 ; end
if ~isfield(Opt,"Accumulate") ; Opt.Accumulate=false ; end
if ~isfield(Opt,"MaskBand") ; Opt.MaskBand="switched" ; end

k=CtrlVar.CurrentRunStepNumber;
N=MUA.Nnodes;

% (10 Oct 2026) Nodes whose constraint state changed in this time step, and in the two previous time steps (see the masking below)
SwitchedNow=setxor(BCs1.hPosNode(:),hPosNodeStartOfStep(:));
[RunInfo,SwitchedRecent,isHistoryValid]=UpdateSwitchHistory(RunInfo,SwitchedNow,N);

dt=F1.dt; dtm1=F0.dtRates; dtm2=Fm1.dtRates;
isPositiveScalar=@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0 ;

Flag=0;
if abs(CtrlVar.theta-0.5)>1e-10
    Flag=1;
elseif ~(isPositiveScalar(dt) && isPositiveScalar(dtm1) && isPositiveScalar(dtm2))
    Flag=2;
elseif numel(F0.h)~=N || numel(F1.h)~=N || numel(F0.dhdt)~=N || numel(Fm1.dhdt)~=N
    Flag=3;
end

if Flag>0
    RunInfo=SetSeries(RunInfo,k,F1.time,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,Flag,NaN,dt);
    if numel(F1.h)==N
        F1=UpdateAccumulated(F1,Opt,[],dt,N);    % only TotalTime is increased
    end
    if Opt.InfoLevel>=1
        Reason=["theta is not 0.5","no second-order predictor available","inconsistent array sizes"];
        fprintf(' Time-discretisation error estimate (h): not calculated (%s). \n',Reason(Flag))
    end
    return
end

%% predictor (unclipped) and corrector
hC=F1.h(:); X=F0.h(:); D=F0.dhdt(:); Dm1=Fm1.dhdt(:);
hP=X+dt*(D+(D-Dm1)*(dtm1+dt)/(dtm2+dtm1));

%% coefficients
tp=dt; t1=-dtm1; t2=-dtm1-dtm2;       % times relative to t(n)
w1=tp*(tp-t2)/(t1*(t1-t2));            % weight of h(t(n-1)) in the extrapolation to t(n+1)
w2=tp*(tp-t1)/(t2*(t2-t1));            % weight of h(t(n-2))
CC=-dt^3/12; c1=-dtm1^3/12; c2=-dtm2^3/12;
CP=dt*(dt+dtm1)*(dt+dtm1+dtm2)/6;
W=w1*c1+w2*(c1+c2);

e=CC*(hC-hP)/(CP-CC-W);

%% masking
Mask=~isfinite(e);
Mask(BCs1.hFixedNode(:))=true;

ActiveNodes=unique([BCs1.hPosNode(:);hPosNodeStartOfStep(:)]);
Mask(ActiveNodes)=true;
if string(Opt.MaskBand)=="active" || ~isHistoryValid
    BandCentre=ActiveNodes;          % previous rule: neighbours of all active nodes
else
    BandCentre=SwitchedRecent;       % neighbours of nodes whose constraint state changed within the last three time steps
end
if ~isempty(BandCentre)
    isEle=any(ismember(MUA.connectivity,BandCentre),2);
    Mask(unique(MUA.connectivity(isEle,:)))=true;
end

hLimit=CtrlVar.ThickMin;
if isfield(CtrlVar,"ThicknessPenalty") && CtrlVar.ThicknessPenalty
    delta=max(CtrlVar.ThickMin,CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.deltaAbs);
    hLimit=CtrlVar.ThickMin+2*delta;
end
Mask=Mask | hC<=hLimit+eps(hLimit) | X<=hLimit+eps(hLimit);

if isprop(F1,"LSFMask") && ~isempty(F1.LSFMask) && isfield(F1.LSFMask,"NodesOut") && numel(F1.LSFMask.NodesOut)==N
    Mask=Mask | F1.LSFMask.NodesOut(:);
end

e(Mask)=NaN;
F1.hTimeDiscretisationErrorEstimate=e;
F1=UpdateAccumulated(F1,Opt,e,dt,N);

%% norms
% (10 Oct 2026) Norms are calculated with the consistent mass matrix M, as sqrt(x'*M*x/A), where x is set to zero at masked nodes and
% A is the area of the region where the estimate is defined (the sum of the rows of M for these nodes). This is the L2 norm over that region,
% divided by the square root of its area, ie an average (rms) value.
M=MUA.M;
if isempty(M) ; M=MassMatrix2D1dof(MUA) ; end
mRow=full(sum(M,2));
v=~isnan(e);
Av=sum(mRow(v));
ValidArea=Av/sum(mRow);

if any(v)
    nrm=@(x) sqrt(max(full(x'*(M*x)),0)/Av);
    e0=zeros(N,1); e0(v)=e(v);
    RMS=nrm(e0);
    iv=find(v); [Max,i]=max(abs(e(v))); iMax=iv(i);
    xMax=MUA.coordinates(iMax,1); yMax=MUA.coordinates(iMax,2);
    Volume=mRow'*e0;          % = integral of e, ie 1'*M*e
    sc=zeros(N,1); sc(v)=e(v)./(Opt.atol+Opt.rtol*abs(hC(v)));
    Scaled=nrm(sc);
    % (10 Oct 2026) Error per unit time, scaled by the tolerances of the accuracy-based time-step limit (see ATSTimeDiscretisationErrorLimit.m)
    atolRate=0.01; rtolRate=0;
    if isfield(CtrlVar,"ATSTimeDiscretisationError")
        if isfield(CtrlVar.ATSTimeDiscretisationError,"atol") ; atolRate=CtrlVar.ATSTimeDiscretisationError.atol ; end
        if isfield(CtrlVar.ATSTimeDiscretisationError,"rtol") ; rtolRate=CtrlVar.ATSTimeDiscretisationError.rtol ; end
    end
    qRate=zeros(N,1); qRate(v)=(e(v)/dt)./(atolRate+rtolRate*abs(hC(v)));
    RateScaled=nrm(qRate);
else
    RMS=NaN; Max=NaN; xMax=NaN; yMax=NaN; Volume=NaN; Scaled=NaN; RateScaled=NaN;
end

RunInfo=SetSeries(RunInfo,k,F1.time,RMS,Max,xMax,yMax,Volume,Scaled,RMS/dt,ValidArea,0,RateScaled,dt);

if Opt.InfoLevel>=1
    fprintf(' Time-discretisation error estimate (h), t=%g (dt=%g): rms=%-9.3g m, max=%-9.3g m at (x,y)=(%g,%g) km, volume=%-9.3g m^3, scaled=%-9.3g, valid area %4.1f%% \n',...
        F1.time,dt,RMS,Max,xMax/1000,yMax/1000,Volume,Scaled,100*ValidArea)
end

end


function RunInfo=SetSeries(RunInfo,k,Time,RMS,Max,xMax,yMax,Volume,Scaled,RMSperUnitTime,ValidArea,Flag,RateScaled,DtStep)

% Stores the values for run step k in the time series in RunInfo.Forward. The series are created, or extended, with NaN if needed (eg for
% RunInfo objects from older restart files).

Names=["hTimeDiscretisationErrorTime","hTimeDiscretisationErrorRMS","hTimeDiscretisationErrorMax","hTimeDiscretisationErrorMaxX",...
    "hTimeDiscretisationErrorMaxY","hTimeDiscretisationErrorVolume","hTimeDiscretisationErrorScaled","hTimeDiscretisationErrorRMSperUnitTime",...
    "hTimeDiscretisationErrorValidArea","hTimeDiscretisationErrorFlag","hTimeDiscretisationErrorRateScaled","hTimeDiscretisationErrorDt"];
Values=[Time,RMS,Max,xMax,yMax,Volume,Scaled,RMSperUnitTime,ValidArea,Flag,RateScaled,DtStep];

for I=1:numel(Names)
    if ~isfield(RunInfo.Forward,Names(I))
        RunInfo.Forward.(Names(I))=NaN(max(k,numel(RunInfo.Forward.time)),1);
    elseif numel(RunInfo.Forward.(Names(I)))<k
        RunInfo.Forward.(Names(I))=[RunInfo.Forward.(Names(I))(:);NaN(k-numel(RunInfo.Forward.(Names(I))),1)];
    end
    RunInfo.Forward.(Names(I))(k)=Values(I);
end

end


function F1=UpdateAccumulated(F1,Opt,e,dt,N)

% (10 Oct 2026) Updates the accumulated time-discretisation error estimate, see the description at the top. The accumulated fields are
% initialised with zeros if they do not exist, or if their size does not agree with the number of nodes. Masked (NaN) values are not
% included. If e is empty (no estimate calculated in this time step), only TotalTime is increased.

if ~Opt.Accumulate
    return
end

A=F1.hTimeDiscretisationErrorAccumulated;
if ~isstruct(A) || ~isfield(A,"Signed") || numel(A.Signed)~=N
    A=struct("Signed",zeros(N,1),"Abs",zeros(N,1),"ValidTime",zeros(N,1),"TotalTime",0);
end

if ~isempty(e)
    v=~isnan(e);
    A.Signed(v)=A.Signed(v)+e(v);
    A.Abs(v)=A.Abs(v)+abs(e(v));
    A.ValidTime(v)=A.ValidTime(v)+dt;
end

A.TotalTime=A.TotalTime+dt;
F1.hTimeDiscretisationErrorAccumulated=A;

end


function [RunInfo,SwitchedRecent,isHistoryValid]=UpdateSwitchHistory(RunInfo,SwitchedNow,N)

% (10 Oct 2026) Keeps the nodes whose constraint state changed in the two previous time steps in RunInfo.Forward, and returns the union of
% these and the nodes that changed in this time step. The history is only valid if it covers the two previous time steps on the same
% mesh (same number of nodes); otherwise it is restarted (eg at the start of a run, or after remeshing).

Name="hTimeDiscretisationErrorSwitchHistory";
if isfield(RunInfo.Forward,Name) && isstruct(RunInfo.Forward.(Name)) && isfield(RunInfo.Forward.(Name),"Nnodes") && RunInfo.Forward.(Name).Nnodes==N
    H=RunInfo.Forward.(Name);
else
    H=struct("Nnodes",N,"Previous",{{}});
end

isHistoryValid=numel(H.Previous)>=2;
SwitchedRecent=unique([SwitchedNow(:);vertcat(H.Previous{:})]);

H.Previous=[{SwitchedNow(:)},H.Previous];
H.Previous=H.Previous(1:min(2,numel(H.Previous)));
RunInfo.Forward.(Name)=H;

end
