
function [UserVar,RunInfo,F1,F0,l,Kuv,Ruv,Lubvb]= uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,Fm1,l)


nargoutchk(8,8)
narginchk(8,8)

%% 1) Calculate uv at the beginning of the time step.  
%
% This is not expected to be
% strictly needed as uv should have been calculated on the basis of the new
% geometry (not the F0 geometry) at the end of last semi-implicit time step
% (Step 4 below). So I could consider getting rid of this, but it should not cost
% too much time.
%
% [UserVar,RunInfo,F0,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,l);

%% 2) Get an explicit estimate of uv

if CtrlVar.InfoLevel>=10 
    fprintf("uvhSemiImplicit:Getting an explicit estimate for u,v and h at t=t1.\n")
end


CtrlVar.StartSemiImplicitWithExtrapolation=false ; % TestIng

if CtrlVar.StartSemiImplicitWithExtrapolation
    
    CtrlVar.ExplicitEstimationMethod="-Adams-Bashforth-";
    [UserVar,RunInfo,F1.ub,F1.vb,F1.ud,F1.vd,F1.h]=ExplicitEstimationForUaFields(UserVar,RunInfo,CtrlVar,MUA,F0,Fm1,BCs,l,BCs,l);
    % The explicit estimate for velocities is the one used for the velocities at the
    % end of the time step when h is calculated implicitly.
    
    % For the velocities at the end of the time step. The explicit estimate for the
    % thickness (h) is the initial guess for h, which is then updated as h is solved
    % implicitly.
    
    CtrlVar.time=CtrlVar.time+CtrlVar.dt;  % I here need the mass balance at the end of the time step, hence must increase t
    
    F1.GF=GL2d(F1.B,F1.S,F1.h,F1.rhow,F1.rho,MUA.connectivity,CtrlVar);
    [UserVar,F1]=GetMassBalance(UserVar,CtrlVar,MUA,F1);
else
    F1=F0;
end

CtrlVar.time=CtrlVar.time-CtrlVar.dt; % and then take it back to t at the beginning.

dadt=(F1.as+F1.ab-(F0.as+F0.ab))/CtrlVar.dt;


if CtrlVar.InfoLevel>=10 
    fprintf("uvhSemiImplicit:Solving for h at t=t1, implicitly.\n")
end

iCount=0 ; Du=inf ; Dv=inf ; Dh=inf ; 
while ( Du>0.1 || Dv>0.1 || Dh> 0.1 ) && iCount < 2 
    
    iCount=iCount+1; 
    h1=F1.h ; u1=F1.ub ; v1=F1.vb;
    
    % 3) calculate new ice thickness implicitly with respect to h.
    [UserVar,RunInfo,F1.h,l]=SSS2dPrognostic(UserVar,RunInfo,CtrlVar,MUA,BCs,l,F0.h,F0.ub,F0.vb,F0.dubdt,F0.dvbdt,F0.as+F0.ab,dadt,F1.ub,F1.vb,F1.as+F1.ab,dadt,F1.dubdt,F1.dvbdt);
    % Make sure to update s and b too.
    [F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);

    % 4) And finally calculate uv based on the new geometry.
    [UserVar,RunInfo,F1,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F1,l);

    du=u1-F1.ub; dv=v1-F1.vb ; dh=h1-F1.h;
    Du=norm(du)/norm(F1.ub) ; Dv=norm(dv)/norm(F1.vb) ; Dh=norm(dh)/norm(F1.h);
    if Du > 0.25 || Dv > 0.25 || Dh> 0.05  % if so, then force a reduction in dt
        RunInfo.Forward.Iterations=max([10*CtrlVar.ATSTargetIterations ; RunInfo.Forward.Iterations]);
        RunInfo.Forward.Iterations=666 ; 
    end
    [Du Dv Dh]
    
end


if CtrlVar.InfoLevel>=10 
    fprintf("uvhSemiImplicit:Solving for uv for t1.\n")
end


F1.dhdt=(F1.h-F0.h)/CtrlVar.dt;
F1.dsdt=(F1.s-F0.s)/CtrlVar.dt;
F1.dbdt=(F1.b-F0.b)/CtrlVar.dt;
F1.dubdt=(F1.ub-F0.ub)/CtrlVar.dt ;
F1.dvbdt=(F1.vb-F0.vb)/CtrlVar.dt;
F1.duddt=(F1.ud-F0.ud)/CtrlVar.dt ;
F1.dvddt=(F1.vd-F0.vd)/CtrlVar.dt;


end
