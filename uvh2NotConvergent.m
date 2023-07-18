
function [UserVar,RunInfo,F1,F0,l0,l1,BCs1,dtOut]=uvh2NotConvergent(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)

dtIn=CtrlVar.dt ;
isF0reset=false;

if CtrlVar.AdaptiveTimeStepping
    dtMin=max(CtrlVar.ATSdtMin,CtrlVar.dtmin)  ;
else
    dtMin=CtrlVar.dtmin ;
end


while CtrlVar.dt > dtMin

    F1.ub=F0.ub ; F1.vb=F0.vb ;F1.ud=F0.ud ;F1.ud=F0.ud ; F1.h=F0.h ; l1=l0 ;
    F1.h(F1.h<CtrlVar.ThickMin)=CtrlVar.ThickMin;

    CtrlVar.dt=CtrlVar.dt/2 ;

    F1.dt=CtrlVar.dt ;  F0.dt=CtrlVar.dt ; dtOut=CtrlVar.dt ;

    fprintf("uvh2NotConvergent: uvh solve did not converge. Reducing time step to dt=%g and try solving again. \n,",CtrlVar.dt)

    [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);


    if RunInfo.Forward.uvhConverged
        fprintf("uvh2NotConvergent: uvh solve converged with dt=%g. \n,",CtrlVar.dt)
        break
    end

    if ~isF0reset  && CtrlVar.dt < dtIn/10
        % OK, I've reduced the original time step by a factor of 10 by now, and still not finding convergence
        %     Will now reset F0
        fprintf("uvh2NotConvergent: Will now reset solution and calulate new uv starting point. \n,")
        F0.ub=F0.ub*0 ; F0.vb=F0.vb*0 ;F0.ud=F0.ud*0 ;F0.ud=F0.ud*0 ;
        F0=StartVelocity(CtrlVar,MUA,BCs1,F0) ;
        [UserVar,RunInfo,F0,l0] = uv(UserVar,RunInfo,CtrlVar,MUA,BCs1,F0,l0);
        isF0reset=true;
    end


end


end