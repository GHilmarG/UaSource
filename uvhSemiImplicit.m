




function [UserVar,RunInfo,F1,F0,l1,Kuv,Ruv,Lubvb,duv1NormVector]= uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) 
                                                             %   uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,Fm1,l) ;


nargoutchk(5,9)
narginchk(8,8)


%% A semi-implicit solver for velocities and thickness, where the thickness is solved implicitly using the theta method, and the velocity explicitly.
%
% Let
%
% $$ F(u,v,h) = 0  $$ 
%
% be the momentum equations, and 
%
%
% $$ M(u,v,h)=0$$ 
%
% the mass conservation equation.
%
% We have the solution for u, v and h at $t=t_0$ and we want to calculate the solution at $t=t_1$ where $t_1=t_0+\Delta t$ 
%
% If we have not done so already, we initially solve for $u_0$ and $v_0$ for a known $h_0$, i.e.
%
% solve $F(u_0,v_0,h_0) = 0$ for $u_0$ and $v_0$ with $h_0$ known.
%
% We then make an initial explicit guess for what the velocities are for $t=t_1$. We call this guess
%
% $u_1^0$ and $v_1^0$, where the upper index is the $u_1$, $v_1$ iterate. 
%
% We then solve $M(u,v,h)=0$ implicitly for $h_1^0$ using $h_0$, $u_0$, $v_0$ and our guesses $u_1^0$ and $v_1^0$ 
%
% Now that we have the estimate $h_1^0$ at the end of the time step, we can use the momentum equation to solve for $u_1$ and $v_1$ using $h_1^0$. 
% 
% This involves solving
%
% $$ F(u_1^1,v_1^1,h_1^0) = 0 $$ 
%
% for $u_1^1$ and  $v_1^1$, which is our new velocity iterate to $t=t_1$.
% 
% We now check by how much the velocity iterate has changed by evaluating
%
%
% $$ e_{uv} := \| (u_1^1,v_1^1) - (u_1^0,v_1^0) \|  $$
%
% If this error is not within a given tolerance, we conclude that our previous estimate, $u_1^0$ amd $v_0^1$ were not good enough, 
% and we recalculate the thickness again implicitly using $u_1^1$ and $v_1^1$.
%
% This results in a new estimate for $h_1$ which we call $h_1^1$ and we evaluate 
%
% $$ e_{h} := \| h_1^1 - h_0 \|  $$
%
%
% We can also check if this is within a given tolerance. 
%
% We then continue to recalculate $u_1^i$ and $v_1^i$ for $i=1,2,3,...$ until a give tolerance is met.
% 
% *The procedure is as follows:* 
%
% For a given initial $h_0$, calculate $u_0$ and $v_0$ using the momentum equation,
% and for 
% 
% $$i=0$$ 
% 
% set 
% 
% $$u_1^i=u_0$$ 
% 
% and
% 
% $$v_0^i=v_0$$ 
% 
% (or use some explicit estimate).
%
% Then do the loop: 
%
% while $e_{uv} > \epsilon$ do: 
%
% $$ i = i+1 $$ 
%
% Find $h_1^i$ by solving the mass-conservation equation implicitly using $u_1^{i-1}$, $v_1^{i-1}$ and $u_0$, $v_0$ and $h_0$
%
% Now find $u_1^{i}$ and $v_1^{i}$ by solving the momentum equation using $h_1^i$. 
%
% Calculate a norm of the changes in the velocity iterate : 
%
% $$ e_{uv} := \| (u_1^i,v_1^i) - (u_1^{i-1},v_1^{i-1}) \|   $$
%
% end 
%
%%




%% If required, calculate uv at the beginning of the time step, ie at t=t0;
if CtrlVar.InitialDiagnosticStep   % if not a restart step, and if not explicitly requested by user, then do not do an initial diagnostic step
    %% diagnostic step, solving for uv.  Always needed at a start of a transient run. Also done if requested by the user.
    CtrlVar.InitialDiagnosticStep=0;

    fprintf(" initial diagnostic step at t=%-.15g \n ",F0.time);

    [UserVar,RunInfo,F0,l1]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs1,F0,l1);
    


end

% [UserVar,F0]=GetCalving(UserVar,CtrlVar,MUA,F0,BCs);  
                                               
%% Get an initial estimate of uv at end of time step, ie at t=t1
CtrlVar.StartSemiImplicitWithExtrapolation=false; % update this later...for the time being I'm not doing an initial explicit estimate for h. 
                                                  % This requires having Fm1 as an input argument. Simple to change at some later stage.


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
    F1=F0; % The starting point for the future state is here simply the current state.
end



%% Now uv at t=t0 is available, as is an initial estimate for uv at t=t1
%
% 1) Solve h1=h(h0,uv0,uv1Est)
% 2) Solve uv1=uv(h1)
% 3) Estimate difference between uv1 and uv1Est. If small, then exit. Otherwise set uv1Est=uv1 and go back to step 1
%
%



duv1NormVector=nan(CtrlVar.uv2h.MaxIterations,1) ;
dh1NormVector=nan(CtrlVar.uv2h.MaxIterations,1) ;

uvItMax=0;
hItMax=0; 
CtrlVar.InfoLevelNonLinIt=0; % here setting the info level to zero, as I give info on uv-h outer iteration and just the number of inner -uv- and -h- iterations

for uv2hIt=1:CtrlVar.uv2h.MaxIterations  % this is the "outer" iteration

    %% 1) calculate new ice thickness implicitly with respect to h: h1=h(h0,uv0,uv1Est)

    if CtrlVar.InfoLevel>=10
        fprintf("uvhSemiImplicit:Solving for changes in h implicitly, going from t=t0=%g to t=t1=%g \n",CtrlVar.time,CtrlVar.time+CtrlVar.dt)
    end


    h1Ahead=F1.h;  % keep this to track changes in the estimate h at t=t1

                                 % now solve for the thickness, this may involve iterating over h, and these are the "inner" h-iterations. 
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
    % Now solve for the velocities. This may involve some uv iterations, and these are the "inner" uv iterations
    [UserVar,RunInfo,F1,l1,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs1,F1,l1);

    uvItMax=max(RunInfo.Forward.uvIterations(CtrlVar.CurrentRunStepNumber),uvItMax);

    du1=F1.ub-ub1Ahead ;
    dv1=F1.vb-vb1Ahead ;
    dh1=F1.h-h1Ahead;

    if ~isempty(F1.LSF)

        % If using the level set, set all differences downstream of calving fronts to zero (this is not "real" ice) 
        if isempty(F1.LSFMask)
            F1.LSFMask=CalcMeshMask(CtrlVar,MUA,F1.LSF,0);
        end

        I=F1.LSFMask.NodesOut ; 
        du1(I)=0 ; dv1(I)=0 ; dh1(I)=0 ;

    end

    duv1Norm=(du1'*MUA.M*du1+dv1'*MUA.M*dv1)/MUA.Area;  % Note: this is actually second power of the norm, not the norm itself
    dh1Norm=dh1'*MUA.M*dh1/MUA.Area;


    duv1NormVector(uv2hIt)=duv1Norm;
    dh1NormVector(uv2hIt)=dh1Norm;
    fprintf("=============     uv-h: Outer iteration %i \t |duv|^2 =%14g \t |dh|^2=%14g    \t inner uv iterations=%2i \t inner h iterations=%2i   \n",...
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

    if duv1Norm< CtrlVar.uv2h.uvTolerance
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
    semilogy(1:CtrlVar.uv2h.MaxIterations,duv1NormVector,DisplayName="$\|\Delta (uv)_1\|$",Marker="o")
    ylabel("$\|\Delta (uv)_1\|$",Interpreter="latex")
    yyaxis right
    semilogy(1:CtrlVar.uv2h.MaxIterations,dh1NormVector,DisplayName="$\|\Delta h_1\|$",Marker="o")
    ylabel("$\|\Delta h_1\|$",Interpreter="latex")
    xlabel("$uv-h$ outer iteration",Interpreter="latex")
    title("$uv-h$ outer iteration residuals",Interpreter="latex")
    subtitle(sprintf("$t$=%g \\quad $\\Delta t$=%g",F1.time,F1.dt),Interpreter="latex")
    legend(Location="best",Interpreter="latex");
end

end
