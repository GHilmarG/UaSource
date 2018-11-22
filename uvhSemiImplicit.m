
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
[UserVar,RunInfo,F0,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,l);

%% 2) Get an explicit estimate of uv

F1=F0;
F1=ExplicitEstimationForUaFields(CtrlVar,F1,F0,Fm1);  
% The explicit estimate for velocities is the one used for the velocities at the
% end of the time step when h is calculated implicitly.

% For the velocities at the end of the time step. The explicit estimate for the
% thickness (h) is the initial guess for h, which is then updated as h is solved
% implicitly.

CtrlVar.time=CtrlVar.time+CtrlVar.dt;  % I here need the mass balance at the end of the time step, hence must increase t

F1.GF=GL2d(F1.B,F1.S,F1.h,F1.rhow,F1.rho,MUA.connectivity,CtrlVar);
[UserVar,F1]=GetMassBalance(UserVar,CtrlVar,MUA,F1);
CtrlVar.time=CtrlVar.time-CtrlVar.dt; % and then take it back to t at the beginning.

dadt=(F1.as+F1.ab-(F0.as+F0.ab))/CtrlVar.dt;

%% 3) calculate new ice thickness implicitly with respect to h.
[F1.h,l]=SSS2dPrognostic(CtrlVar,MUA,BCs,l,F0.h,F0.ub,F0.vb,F0.dubdt,F0.dvbdt,F0.as+F0.ab,dadt,F1.ub,F1.vb,F1.as+F1.ab,dadt,F1.dubdt,F1.dvbdt);
% Make sure to update s and b too.
[F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);
%[F1.b,F1.s,F1.h]=Calc_bs_From_hBS(F1.h,F1.S,F1.B,F1.rho,F1.rhow,CtrlVar,MUA.coordinates);

%% 4) And finally calculate uv based on the new geometry.

[UserVar,RunInfo,F1,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F1,l);

% At the end of the time step, the velocities are fully consistent with he
% thickness at the end of the time step.
%
% The only inconsitency is that the velocities at the end of the time step used
% to calculate the new thickness are those from the explicit estimate.

F1.dhdt=(F1.h-F0.h)/CtrlVar.dt;
F1.dsdt=(F1.s-F0.s)/CtrlVar.dt;
F1.dbdt=(F1.b-F0.b)/CtrlVar.dt;
F1.dubdt=(F1.ub-F0.ub)/CtrlVar.dt ;
F1.dvbdt=(F1.vb-F0.vb)/CtrlVar.dt;
F1.duddt=(F1.ud-F0.ud)/CtrlVar.dt ;
F1.dvddt=(F1.vd-F0.vd)/CtrlVar.dt;


end
