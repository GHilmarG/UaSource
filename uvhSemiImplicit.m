




function [UserVar,RunInfo,F1,F0,l,Kuv,Ruv,Lubvb,duv1NormVector]= uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,Fm1,l)


nargoutchk(9,9)
narginchk(8,8)

%% If required, calculate uv at the beginning of the time step, ie at t=t0;
if CtrlVar.InitialDiagnosticStep   % if not a restart step, and if not explicitly requested by user, then do not do an initial diagnostic step
    %% diagnostic step, solving for uv.  Always needed at a start of a transient run. Also done if requested by the user.
    CtrlVar.InitialDiagnosticStep=0;

    fprintf(" initial diagnostic step at t=%-.15g \n ",CtrlVar.time);

    [UserVar,RunInfo,F0,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,l);
    


end

% [UserVar,F0]=GetCalving(UserVar,CtrlVar,MUA,F0,BCs);  
                                               
%% Get an initial estimate of uv at end of time step, ie at t=t1
CtrlVar.StartSemiImplicitWithExtrapolation=false; % update this later

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


%% Now uv at t=t0 is available, as is an initial estimate for uv at t=t1
% 
% 1) Solve h1=h(h0,uv0,uv1Est)
% 2) Solve uv1=uv(h1)
% 3) Estimate difference between uv1 and uv1Est. If small, then exit. Otherwise set uv1Est=uv1 and go back to step 1
%
%
uvItMax=5; Tolerance= 1e-5 ; 
duv1NormVector=nan(uvItMax,1) ;

for uvIt=1:uvItMax

    %% calculate new ice thickness implicitly with respect to h: h1=h(h0,uv0,uv1Est)

    h1Ahead=F1.h; 
    [UserVar,RunInfo,F1.h,l]=SSS2dPrognostic(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);

    % Make sure to update s and b too.
    [F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);

    %% Calculate uv for the new thickness: uv1=uv(h1)
    

    ub1Ahead=F1.ub;  vb1Ahead=F1.vb;
    [UserVar,RunInfo,F1,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F1,l);


    du1=F1.ub-ub1Ahead ;
    dv1=F1.vb-vb1Ahead ;


    duv1Norm=norm([du1;dv1])/sqrt(2*MUA.Nnodes) ;
    dh1Norm=norm(F1.h-h1Ahead)/sqrt(MUA.Nnodes) ;
    duv1NormVector(uvIt)=duv1Norm; 

    FT=sprintf("duv1-%i",uvIt) ;
    UaPlots(CtrlVar,MUA,F1,[du1 dv1],FigureTitle=FT)
    title(FT)
    subtitle(sprintf("t=%g   dt=%g",CtrlVar.time,CtrlVar.dt),Interpreter="latex")

    fprintf("duvh1Norm=%g \t dh1Norm=%g \n",duv1Norm,dh1Norm)
    if duv1Norm< Tolerance
        break
    end


end

if CtrlVar.InfoLevel>=10
    fprintf("uvhSemiImplicit: Now solved for uv and h at t=t1.\n")
end



end
