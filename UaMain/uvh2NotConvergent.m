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
% Called when a uvh solve did not converge. Two rescue branches are tried, in this order:
%
%   1) (12 Nov 2026) Increase the transition length l of the grounding-line mask (F0.lGL, F1.lGL), see TransitionLengthRescue below.
%      Only tried if the failure can be related to l, ie if integration points lie within the transition of the mask, and only if
%      CtrlVar.GroundingLineTransition.Rescue.isOn. l is doubled (at most to lmax) at the nodes of those elements, and the solve repeated
%      with the same time step and the same active set, until it converges. The increased l is kept for the rest of this time step (it is returned in F0), and
%      at the start of the next time step l is again calculated from GroundingLineTransitionLength.m. If
%      CtrlVar.GroundingLineTransition.Rescue.ReturnToDefault, a final solve with the original l is done, starting from the solution
%      obtained with the increased l, and that solution is accepted if it converges.
%
%   2) Repeatedly halve the time step and solve again until either the solve converges, or the minimum allowed time step is reached.
%      This is only done if branch 1 was not tried or did not converge.
%
% RunInfo.Forward.uvhTransitionLengthRescued is true on return if the solve converged in branch 1 (ie with the time step unchanged).
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

%% (12 Nov 2026) Rescue branch 1: increase the transition length of the grounding-line mask
[UserVar,RunInfo,F1,F0,l1,BCs1,Rescued]=TransitionLengthRescue(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1) ;
RunInfo.Forward.uvhTransitionLengthRescued=Rescued ;
if Rescued
    return
end

%% Rescue branch 2: reduce the time step

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [UserVar,RunInfo,F1,F0,l1,BCs1,Rescued]=TransitionLengthRescue(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)
%%
% (12 Nov 2026) Rescue branch 1 of uvh2NotConvergent: increase the transition length l of the grounding-line mask.
%
% Why: when the grounding line moves such that a group of integration points passes through the transition of the mask, the
% feedback between grounding (basal drag), velocity and thickness at these integration points can fold the discrete problem:
% within the transition the mask changes by order one over a thickness change of order l, and if l is small compared with the
% change of h-hf at these points caused by their own grounding, the residual is not monotonic in h-hf there. The Newton
% direction then points away from the root on one side of the fold, and the iteration ends at the fold, where the Jacobian is
% (nearly) singular in that direction and the line search fails.
% 
%
% The rescue is only tried if integration points lie within the transition, ie if |h-hf-Hh0| < k*l at an integration point at
% the start of the time step (F0) or at the last iterate (F1), with k=CtrlVar.GroundingLineTransition.Rescue.k. Otherwise the
% mask is effectively a step function at all integration points, l plays no role in the failure, and the time step is reduced
% instead (rescue branch 2).
%
% l is doubled, at most to lmax, at the nodes of the elements with integration points within the transition (the set is
% extended after each failed attempt), and the solve is repeated from the start of the time step, with the same time step and
% with the active set of the failed solve (BCs1 on input): unlike in branch 2, where the active set is reset to the start of
% the time step, the time step is unchanged and that active set is therefore still appropriate.
% 
% At most CtrlVar.GroundingLineTransition.Rescue.MaxIncreases attempts are made. If the solve converges, and if
% CtrlVar.GroundingLineTransition.Rescue.ReturnToDefault is true, a further solve is done with the original l, starting from
% that solution. If this converges, the solution with the original l is accepted, ie the rescue only changes how the solution
% is found; otherwise the solution with the increased l is kept.
%
% On return F0.lGL and F1.lGL are the transition lengths of the accepted solution, which are then used for the rest of this
% time step (the caller keeps F0). At the start of the next time step l is again calculated from
% GroundingLineTransitionLength.m. If the rescue is not tried or does not converge, F1, l1, BCs1 and the transition lengths
% are returned as on input.
%
% Each rescue is recorded in RunInfo.Forward.uvhTransitionLengthRescueLog, one row per rescue with the columns: run step,
% time, dt, number of increases, largest l, number of nodes where l was increased, converged (0/1), returned to the original
% l (0/1).
%%

Rescued=false ;
GLT=CtrlVar.GroundingLineTransition ;
if ~isfield(GLT,"Rescue") || ~GLT.Rescue.isOn
    return
end
Opt=GLT.Rescue ;

lInput=F0.lGL ;   % transition length on input ([] means the scalar lmin)
lInputNod=lInput ; if numel(lInputNod)~=MUA.Nnodes ; lInputNod=GLT.lmin+zeros(MUA.Nnodes,1) ; end
lNod=lInputNod ;

Ele=ElementsWithinTransition(CtrlVar,MUA,F0,lNod,Opt.k) | ElementsWithinTransition(CtrlVar,MUA,F1,lNod,Opt.k) ;
if ~any(Ele)
    fprintf("uvh2NotConvergent: no integration points within the grounding-line transition (|h-hf|<%g l). The transition length is not increased. \n",Opt.k)
    return
end

F1Input=F1 ; l1Input=l1 ; BCs1Input=BCs1 ; RunInfoForwardInput=RunInfo.Forward ;
nIncreases=0 ;

while nIncreases < Opt.MaxIncreases

    Nodes=unique(MUA.connectivity(Ele,:)) ;
    lNew=lNod ; lNew(Nodes)=min(GLT.lmax,2*lNod(Nodes)) ;
    if all(lNew==lNod)
        fprintf("uvh2NotConvergent: the transition length is already at lmax=%g at all %i nodes concerned. \n",GLT.lmax,numel(Nodes))
        break
    end
    lNod=lNew ; nIncreases=nIncreases+1 ;

    % start again from the start of the time step, with the active set of the failed solve
    BCs1=BCs1Input ;
    F1.ub=F0.ub ; F1.vb=F0.vb ; F1.ud=F0.ud ; F1.vd=F0.vd ; F1.h=F0.h ; l1=l0 ;
    F1.h(F1.h<CtrlVar.ThickMin)=CtrlVar.ThickMin ;
    F1.h(BCs1.hPosNode)=CtrlVar.ThickMin ;
    F0.lGL=lNod ; F1.lGL=lNod ;

    fprintf("uvh2NotConvergent: uvh solve did not converge. Increasing the grounding-line transition length at %i nodes (to at most %g m, attempt %i) and solving again with dt=%g. \n",...
        nnz(lNod>lInputNod),max(lNod),nIncreases,F1.dt)

    [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;

    if RunInfo.Forward.uvhConverged
        Rescued=true ;
        break
    end

    Ele=Ele | ElementsWithinTransition(CtrlVar,MUA,F1,lNod,Opt.k) ;   % include integration points that have entered the (now wider) transition

end

ReturnedToDefault=false ;
if Rescued

    fprintf("uvh2NotConvergent: uvh solve converged with the increased transition length (at most %g m). \n",max(lNod))

    if Opt.ReturnToDefault
        F1r=F1 ; l1r=l1 ; BCs1r=BCs1 ; RunInfoForwardr=RunInfo.Forward ;
        F0.lGL=lInput ; F1.lGL=lInput ;
        fprintf("uvh2NotConvergent: solving again with the original transition length, starting from this solution. \n")
        [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;
        if RunInfo.Forward.uvhConverged
            ReturnedToDefault=true ;
            fprintf("uvh2NotConvergent: uvh solve converged with the original transition length. \n")
        else
            fprintf("uvh2NotConvergent: uvh solve did not converge with the original transition length. The solution with the increased transition length is used. \n")
            F1=F1r ; l1=l1r ; BCs1=BCs1r ; F0.lGL=lNod ; F1.lGL=lNod ;
            RunInfo.Forward=RunInfoForwardr ;
        end
    end

else

    fprintf("uvh2NotConvergent: uvh solve did not converge with the increased transition length either. \n")
    F1=F1Input ; l1=l1Input ; BCs1=BCs1Input ; F0.lGL=lInput ; F1.lGL=lInput ;
    RunInfo.Forward=RunInfoForwardInput ;

end

Row=[CtrlVar.CurrentRunStepNumber CtrlVar.time F1.dt nIncreases max(lNod) nnz(lNod>lInputNod) Rescued ReturnedToDefault] ;
if ~isfield(RunInfo.Forward,"uvhTransitionLengthRescueLog") || isempty(RunInfo.Forward.uvhTransitionLengthRescueLog)
    RunInfo.Forward.uvhTransitionLengthRescueLog=Row ;
else
    RunInfo.Forward.uvhTransitionLengthRescueLog(end+1,:)=Row ;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ele=ElementsWithinTransition(CtrlVar,MUA,F,lNod,k)
% (12 Nov 2026) Elements with at least one integration point within the transition of the grounding-line mask, |h-hf-Hh0| < k*l, with
% h-hf and l interpolated from the nodes to the integration points (as the transition length in IntegrationPointkH.m).
hf=F.rhow.*(F.S-F.B)./F.rho ;
dnod=reshape(F.h(MUA.connectivity)-hf(MUA.connectivity),MUA.Nele,MUA.nod) ;
lnod=reshape(lNod(MUA.connectivity),MUA.Nele,MUA.nod) ;
Ele=false(MUA.Nele,1) ;
for Iint=1:MUA.nip
    fun=shape_fun(Iint,2,MUA.nod,MUA.points) ;
    Ele=Ele | abs(dnod*fun-CtrlVar.Hh0) < k*(lnod*fun) ;
end
end
