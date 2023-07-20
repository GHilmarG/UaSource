
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



if CtrlVar.InitialDiagnosticStep   % if not a restart step, and if not explicitly requested by user, then do not do an inital dignostic step
    %% diagnostic step, solving for uv.  Always needed at a start of a transient run. Also done if asked by the user.
    CtrlVar.InitialDiagnosticStep=0;

    fprintf(CtrlVar.fidlog,' initial diagnostic step at t=%-.15g \n ',CtrlVar.time);

    [UserVar,RunInfo,F0,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,l);
    


end

[UserVar,F0]=GetCalving(UserVar,CtrlVar,MUA,F0,BCs);  % Level Set  
                                                      % This is the level
                                                      % set at the
                                                      % beginning of the
                                                      % time increment.
                                                      % Currently the level
                                                      %-set is not solved
                                                      % implicitly together
                                                      % with uv. 
                       
                                               


CtrlVar.StartSemiImplicitWithExtrapolation=false ; % TestIng: Hm, for some reason this extrapolation appears to 
                                                  % produce some big
                                                  % oscillations in some
                                                  % runs of, for example,
                                                  % Thwaites, at least in
                                                  % the beginnings of those
                                                  % runs. 
                                                  %
if CtrlVar.StartSemiImplicitWithExtrapolation

    CtrlVar.ExplicitEstimationMethod="-Adams-Bashforth-";
    
    F1=F0;
    [UserVar,RunInfo,F1.ub,F1.vb,F1.ud,F1.vd,F1.h]=ExplicitEstimationForUaFields(UserVar,RunInfo,CtrlVar,MUA,F0,Fm1,BCs,l,BCs,l);
    % The explicit estimate for velocities is the one used for the velocities at the
    % end of the time step when h is calculated implicitly.

    % For the velocities at the end of the time step. The explicit estimate for the
    % thickness (h) is the initial guess for h, which is then updated as h is solved
    % implicitly.

else
    F1=F0;
end


if CtrlVar.InfoLevel>=10
    fprintf("uvhSemiImplicit:Solving for h at t=t1, implicitly.\n")
end




%% calculate new ice thickness implicitly with respect to h.

[UserVar,RunInfo,F1.h,l]=SSS2dPrognostic(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);

% Make sure to update s and b too.
[F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);

%% Calculate uv for the new thickness
%
[UserVar,RunInfo,F1,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F1,l);




%  du=u1-F1.ub; dv=v1-F1.vb ; dh=h1-F1.h;
%  Du=norm(du)/norm(F1.ub) ; Dv=norm(dv)/norm(F1.vb) ; Dh=norm(dh)/norm(F1.h);
%  if Du > 0.25 || Dv > 0.25 || Dh> 0.05  % if so, then force a reduction in dt
%      %RunInfo.Forward.uvhIterations=max([10*CtrlVar.ATSTargetIterations ; RunInfo.Forward.Iterations]);
%      RunInfo.Forward.uvhIterations(CtrlVar.CurrentRunStepNumber)=666;
%  end
% [Du Dv Dh]

% end


if CtrlVar.InfoLevel>=10
    fprintf("uvhSemiImplicit:Solving for uv for t1.\n")
end



end
