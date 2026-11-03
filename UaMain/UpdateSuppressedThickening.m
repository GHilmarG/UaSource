function F=UpdateSuppressedThickening(CtrlVar,MUA,F,BCs,RunInfo)

%% Accumulated suppressed thickening of thickness-constrained nodes (cumulative release criterion of the active-set method)
%
%   F=UpdateSuppressedThickening(CtrlVar,MUA,F,BCs,RunInfo)
%
% (31 Oct 2026) Called once per accepted time step (Ua2D.m), with the converged solution F and the final active set BCs.hPosNode of that
% time step. For each node in the active set, F.hSuppressedThickening (S, units of distance) is the net thickening that the constraint
% has suppressed since the node became constrained:
%
%   S <- max(0, S - ah*dt)        (ah: equivalent mass-balance rate of the constraint, negative where it removes ice,
%                                  see ThicknessConstraintReactions.m)
%
% capped at the release threshold (see below), and S=0 at all other nodes, including nodes released at the end of the time step. A constrained node is released once
% max(0,S-ah*dt) > CtrlVar.ActiveSet.ReleaseAlpha*ThickMin (see ActiveSetUpdate.m, CtrlVar.ActiveSet.ReleaseRule="cumulative"), ie once
% the ice there would have grown by ReleaseAlpha*ThickMin. As S is a time integral of the rate, this does not depend on the time step:
% real accumulation, however small, always leads to the release of the node eventually, while contributions of alternating sign
% (round-off) cancel, and periods of thinning reset S.
%
% The rates ah are those of the last active-set update of the time step (RunInfo.Forward.ActiveSetah, set in ActiveSetUpdate.m from the
% active set and multipliers of the final solve). They are not recomputed here from BCs and l, as these need not match at the end of the
% time step (deferred releases remove nodes from the active set without a further solve).
%
%%

N=MUA.Nnodes ;
S=F.hSuppressedThickening ;
if numel(S)~=N ; S=zeros(N,1) ; end   % eg not yet defined, or from a restart file written before 31 Oct 2026
Snew=zeros(N,1) ;

ah=[] ;
if isfield(RunInfo.Forward,"ActiveSetah") ; ah=RunInfo.Forward.ActiveSetah ; end
% S is capped at the release threshold alpha*hScale: once S exceeds it, the node is released (if the release then fails, eg because
% the node is locked, it stays constrained), so larger values carry no further information. Without the cap, S could grow without
% bound at nodes that remain constrained despite a thickening tendency (eg in the ice-free area behind a calving front, where ice flows
% in and the release fails each time step), and these nodes would then still be released after the tendency has reversed.
alpha=1e-3 ; if isfield(CtrlVar,"ActiveSet") && isfield(CtrlVar.ActiveSet,"ReleaseAlpha") ; alpha=CtrlVar.ActiveSet.ReleaseAlpha ; end
hScaleRelease=CtrlVar.ThickMin ;   % as in ActiveSetUpdate.m
if isfield(CtrlVar,"ThicknessPenaltyMassBalanceFeedbackSoftPlus") && isfield(CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus,"deltaAbs")
    hScaleRelease=max(hScaleRelease,CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.deltaAbs) ;
end

if CtrlVar.ThicknessConstraints && ~isempty(BCs.hPosNode) && numel(ah)==N
    I=BCs.hPosNode(:) ;
    Snew(I)=min(max(0,S(I)-ah(I)*F.dt),alpha*hScaleRelease) ;
end

F.hSuppressedThickening=Snew ;

end
