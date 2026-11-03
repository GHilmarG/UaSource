
function [UserVar,RunInfo,F1,F0,l0,l1,BCs1,dtOut]=uvh2NotConvergent(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1,BCsStartOfStep)

% (2 Nov 2026) BCsStartOfStep (optional): the boundary conditions, including the active set of thickness constraints, at the start of the
% time step. If given, each new attempt with a reduced time step starts from this active set, and not from the active set left by the failed
% attempt. In a failed attempt with a large time step, many nodes can become (deeply) negative and be added to the active set; keeping these
% constraints for the smaller time step held nodes at the minimum thickness that would otherwise keep their ice (Thule: up to 2300 m of ice
% lost at individual nodes), and the cumulative release criterion releases them only gradually.
if nargin<10
    BCsStartOfStep=[] ;
end

%%
%
% Called when a uvh solve did not converge. Repeatedly halves the time step and solves again until either the solve converges, or
% the minimum allowed time step is reached.
%
% The time step is F1.dt, and all decisions here are based on F1.dt. CtrlVar.dt is only a copy of F1.dt, kept for compatibility with
% the rest of the code. On return F1.dt, F0.dt and dtOut are identical and equal to the last time step used.
%
% It is the responsibility of the caller to continue with this time step, i.e. to also set CtrlVar.dt=dtOut. Otherwise subsequent
% solves will be done with the original time step.
%
%%

if ~isfield(CtrlVar,"NeverChangePrescribedTimeStep")
    CtrlVar.NeverChangePrescribedTimeStep=false;
end

dtIn=F1.dt ;
dtOut=dtIn ;
CtrlVar.dt=F1.dt ;     % F1.dt is the master

if CtrlVar.NeverChangePrescribedTimeStep

    fprintf("uvh2NotConvergent: The uvh solve did not converge. But dt is not allowed to be reduced because user has set CtrlVar.NeverChangePrescribedTimeStep to true. \n")
    return

end

isF0reset=false;

if CtrlVar.AdaptiveTimeStepping
    dtMin=max(CtrlVar.ATSdtMin,CtrlVar.dtmin)  ;
else
    dtMin=CtrlVar.dtmin ;
end


if  F1.dt <= dtMin

    fprintf("uvh2NotConvergent: The uvh solve did not converge. But dt can not be reduced any further as it is already set to the minimum allowed time step of %g. \n",dtMin)

else


    while F1.dt > dtMin

        if ~isempty(BCsStartOfStep)   % (2 Nov 2026) start again from the active set at the start of the time step
            nFailed=numel(BCs1.hPosNode) ; BCs1=BCsStartOfStep ;
            RunInfo.Forward.ActiveSetLockedNodes=[] ; RunInfo.Forward.ActiveSetNoFurtherReleases=false ; RunInfo.Forward.ActiveSetah=[] ;
            fprintf("uvh2NotConvergent: active set reset to the start of the time step (%i thickness constraints, against %i after the failed attempt). \n",numel(BCs1.hPosNode),nFailed)
        end
        F1.ub=F0.ub ; F1.vb=F0.vb ; F1.ud=F0.ud ; F1.vd=F0.vd ; F1.h=F0.h ; l1=l0 ;
        F1.h(F1.h<CtrlVar.ThickMin)=CtrlVar.ThickMin;
        F1.h(BCs1.hPosNode)=CtrlVar.ThickMin;       % make sure that the starting point is feasible with respect to the active set

        F1.dt=max(F1.dt/2,dtMin) ;     % never go below the minimum time step
        F0.dt=F1.dt ;
        CtrlVar.dt=F1.dt ;             % copy of F1.dt, kept for compatibility
        F0.lGL=GroundingLineTransitionLength(CtrlVar,MUA,F0) ; F1.lGL=F0.lGL ;   % (2 Nov 2026) transition length for the reduced time step
        dtOut=F1.dt ;

        fprintf("uvh2NotConvergent: uvh solve did not converge. Reducing time step to dt=%g and try solving again. \n",F1.dt)

        [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);


        if RunInfo.Forward.uvhConverged
            fprintf("uvh2NotConvergent: uvh solve converged with dt=%g. \n",F1.dt)
            break
        end

        if ~isF0reset  && F1.dt < dtIn/10
            % OK, I've reduced the original time step by a factor of 10 by now, and still not finding convergence
            %     Will now reset F0
            fprintf("uvh2NotConvergent: Will now reset solution and calculate new uv starting point. \n")
            F0.ub=F0.ub*0 ; F0.vb=F0.vb*0 ; F0.ud=F0.ud*0 ; F0.vd=F0.vd*0 ;
            F0=StartVelocity(CtrlVar,MUA,BCs1,F0) ;
            [UserVar,RunInfo,F0,l0] = uv(UserVar,RunInfo,CtrlVar,MUA,BCs1,F0,l0);
            isF0reset=true;
        end
    end

    if ~RunInfo.Forward.uvhConverged
        fprintf("uvh2NotConvergent: uvh solve did not converge, even with dt=%g. \n",F1.dt)
    end

end

end
