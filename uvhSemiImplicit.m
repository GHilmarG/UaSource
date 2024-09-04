




function [UserVar,RunInfo,F1,F0,l1,Kuv,Ruv,Lubvb,duv1NormVector]= uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) 
                                                             %   uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,Fm1,l) ;


nargoutchk(9,9)
narginchk(8,8)

%% If required, calculate uv at the beginning of the time step, ie at t=t0;
if CtrlVar.InitialDiagnosticStep   % if not a restart step, and if not explicitly requested by user, then do not do an initial diagnostic step
    %% diagnostic step, solving for uv.  Always needed at a start of a transient run. Also done if requested by the user.
    CtrlVar.InitialDiagnosticStep=0;

    fprintf(" initial diagnostic step at t=%-.15g \n ",F0.time);

    [UserVar,RunInfo,F0,l1]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs1,F0,l1);
    


end

% [UserVar,F0]=GetCalving(UserVar,CtrlVar,MUA,F0,BCs);  
                                               
%% Get an initial estimate of uv at end of time step, ie at t=t1
% CtrlVar.StartSemiImplicitWithExtrapolation=false; % update this later
% 
% if CtrlVar.StartSemiImplicitWithExtrapolation
% 
%     CtrlVar.ExplicitEstimationMethod="-Adams-Bashforth-";
% 
%     F1=F0;
%     [UserVar,RunInfo,F1.ub,F1.vb,F1.ud,F1.vd,F1.h]=ExplicitEstimationForUaFields(UserVar,RunInfo,CtrlVar,MUA,F0,Fm1,BCs,l,BCs,l);
%     % The explicit estimate for velocities is the one used for the velocities at the
%     % end of the time step when h is calculated implicitly.
% 
%     % For the velocities at the end of the time step. The explicit estimate for the
%     % thickness (h) is the initial guess for h, which is then updated as h is solved
%     % implicitly.
% 
% else
%     F1=F0;
% end



%% Now uv at t=t0 is available, as is an initial estimate for uv at t=t1
%
% 1) Solve h1=h(h0,uv0,uv1Est)
% 2) Solve uv1=uv(h1)
% 3) Estimate difference between uv1 and uv1Est. If small, then exit. Otherwise set uv1Est=uv1 and go back to step 1
%
%
uv2hItMax=15; Tolerance= 1e-5 ;
duv1NormVector=nan(uv2hItMax,1) ;
dh1NormVector=nan(uv2hItMax,1) ;

uvItMax=0;
hItMax=0; 
for uv2hIt=1:uv2hItMax

    %% 1) calculate new ice thickness implicitly with respect to h: h1=h(h0,uv0,uv1Est)

    if CtrlVar.InfoLevel>=10
        fprintf("uvhSemiImplicit:Solving for changes in h implicitly, going from t=t0=%g to t=t1=%g \n",CtrlVar.time,CtrlVar.time+CtrlVar.dt)
    end


    h1Ahead=F1.h;
    

   [UserVar,RunInfo,h1,l1,BCs1]=MassContinuityEquationNewtonRaphsonThicknessContraints(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ; 
    F1.h=h1; 

    hItMax=max(RunInfo.Forward.hIterations(CtrlVar.CurrentRunStepNumber),hItMax); 

    %% 2)  Update s,h and GF
    [F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);

    %% 3) Calculate uv for the new thickness: uv1=uv(h1)

    if CtrlVar.InfoLevel>=10
        fprintf("uvhSemiImplicit:Solving for uv at t=t1=%g.\n",CtrlVar.time+CtrlVar.dt)
    end


    ub1Ahead=F1.ub;  vb1Ahead=F1.vb;
    [UserVar,RunInfo,F1,l1,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs1,F1,l1);

    uvItMax=max(RunInfo.Forward.uvIterations(CtrlVar.CurrentRunStepNumber),uvItMax);

    du1=F1.ub-ub1Ahead ;
    dv1=F1.vb-vb1Ahead ;
    dh1=F1.h-h1Ahead; 

    if ~isempty(F1.LSF)

        if isempty(F1.LSFMask)
            F1.LSFMask=CalcMeshMask(CtrlVar,MUA,F1.LSF,0);
        end

        I=F1.LSFMask.NodesIn ;  Nnodes=numel(find(I));

        duv1Norm=norm([du1(I);dv1(I)])/sqrt(2*Nnodes) ;
        dh1Norm=norm(dh1)/sqrt(Nnodes) ; 

    else

        duv1Norm=norm([du1;dv1])/sqrt(2*MUA.Nnodes) ;
        dh1Norm=norm(dh1)/sqrt(MUA.Nnodes) ;

    end



    duv1NormVector(uv2hIt)=duv1Norm;
    dh1NormVector(uv2hIt)=dh1Norm;
    fprintf("\n =============     uv-h: Outer iterations %i \t |duv| =%g \t |dh|=%g    \t uv-iterations=%i \t h-iterations=%i  =============  \n",...
        uv2hIt,duv1Norm,dh1Norm,RunInfo.Forward.uvIterations(CtrlVar.CurrentRunStepNumber),RunInfo.Forward.hIterations(CtrlVar.CurrentRunStepNumber))


    if CtrlVar.InfoLevelNonLinIt>= 5 && CtrlVar.doplots

        Tflow=FindOrCreateFigure("flow"); clf(Tflow)
        Ti=tiledlayout('flow');
        nexttile
        UaPlots(CtrlVar,MUA,F1,[du1 dv1],GetRidOfValuesDownStreamOfCalvingFronts=true,CreateNewFigure=false) ;
        title(sprintf("Change in velocites"),Interpreter="latex")
        %subtitle(sprintf("inner iterations: \\#$h$=%i  \\#$uv$=%i",...
        %    RunInfo.Forward.hIterations(CtrlVar.CurrentRunStepNumber),RunInfo.Forward.uvIterations(CtrlVar.CurrentRunStepNumber)),...
        %    Interpreter="latex")

        nexttile
        cbar=UaPlots(CtrlVar,MUA,F1,dh1,GetRidOfValuesDownStreamOfCalvingFronts=true,CreateNewFigure=false) ;
        title(sprintf("Change in ice thickness"),Interpreter="latex")
        title(cbar,"(m)",Interpreter="latex")
        %subtitle(sprintf("inner iterations: \\#$h$=%i  \\#$uv$=%i",...
        %    RunInfo.Forward.hIterations(CtrlVar.CurrentRunStepNumber),RunInfo.Forward.uvIterations(CtrlVar.CurrentRunStepNumber)),...
        %    Interpreter="latex")

        Ti.Title.Interpreter="latex";
        Ti.Title.String=sprintf("$uv-h$ outer iteration %i at $t$=%g and $\\Delta t$=%g : total inner iterations: $h$=%i,  $uv$=%i",...
            uv2hIt,F1.time,F1.dt,RunInfo.Forward.hIterations(CtrlVar.CurrentRunStepNumber),RunInfo.Forward.uvIterations(CtrlVar.CurrentRunStepNumber)) ;
        Ti.Title.Color="blue";

    end

    if duv1Norm< Tolerance
        break
    end

 

end

% Here I return the max number of uv and h iterations during the outer uv-h iteration
RunInfo.Forward.uvIterations(CtrlVar.CurrentRunStepNumber)=uvItMax;
RunInfo.Forward.hIterations(CtrlVar.CurrentRunStepNumber)=hItMax;
RunInfo.Forward.uv2hIterations(CtrlVar.CurrentRunStepNumber)=uv2hIt;

if CtrlVar.InfoLevel>=10

    fprintf("uvhSemiImplicit: Now solved for uv and h at t=t1.\n")

end

if CtrlVar.InfoLevelNonLinIt>= 5 && CtrlVar.doplots
    FindOrCreateFigure("uv-h outer iteration")
    yyaxis left
    semilogy(1:uv2hItMax,duv1NormVector,DisplayName="$\|\Delta (uv)_1\|$",Marker="o")
    ylabel("$\|\Delta (uv)_1\|$",Interpreter="latex")
    yyaxis right
    semilogy(1:uv2hItMax,dh1NormVector,DisplayName="$\|\Delta h_1\|$",Marker="o")
    ylabel("$\|\Delta h_1\|$",Interpreter="latex")
    xlabel("$uv-h$ iteration",Interpreter="latex")
    title("$uv-h$ outer iteration residuals",Interpreter="latex")
    subtitle(sprintf("t=%g dt=%g",F1.time,F1.dt),Interpreter="latex")
    legend(Location="best",Interpreter="latex");
end

end
