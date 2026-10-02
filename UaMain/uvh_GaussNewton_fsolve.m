function [UserVar,RunInfo,F1,l1,BCs1]=uvh_GaussNewton_fsolve(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)

%%
%   [UserVar,RunInfo,F1,l1,BCs1]=uvh_GaussNewton_fsolve(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)
%
% Fully-implicit solution of the uvh problem (SSTREAM) using MATLAB's fsolve (Optimization Toolbox).
%
% This is an alternative to SSTREAM_TransientImplicit and has the same call as uvh2D. It is intended to replace the call to
% SSTREAM_TransientImplicit within uvh2D, so that everything in uvhRootFinding (and uvh2NotConvergent) that surrounds the
% call to uvh2D is retained.
%
% The uvh problem is a root-finding problem with linear constraints:
%
%   R(x)=0 subject to Aeq x = beq,    x=[u;v;h]
%
% where R is the finite-element residual as returned by uvhAssembly and K=dR/dx is its Jacobian. The constraints (boundary
% conditions and ties) are eliminated by writing
%
%   x = xp + N y
%
% where xp is a particular solution and the columns of N are an orthonormal basis of the null space of Aeq. This is done
% by PairwiseConstraintsNullSpace.m, which is possible because the rows of Aeq only contain fixed values and pairwise ties.
%
% The problem that fsolve then solves is the square and unconstrained system
%
%   f(y) = scale * N' R(xp+N y) = 0,        with Jacobian J = scale * N' K N
%
% Note that N'R is the residual projected onto the space of admissible variations, i.e. the part of R that is not balanced
% by the constraints (the reaction forces are not part of it). A root of f is therefore a solution to the constrained root
% problem, and it is the same problem as solved in SSTREAM_TransientImplicit. (Minimising norm(R) directly subject to the
% constraints would be a different problem, because the rows of R that correspond to constrained dofs are not zero.)
%
% With the default algorithm, trust-region-dogleg, fsolve uses the Newton step if it lies within the trust region, and
% otherwise a combination of the Newton and the steepest-descent direction of 0.5*norm(f)^2.
%
% Each time the residual is evaluated, F1.b and F1.s are first recalculated from F1.h (as is done within the Newton loop in
% SSTREAM_TransientImplicit), and if CtrlVar.MassBalanceGeometryFeedback>0 the surface and basal mass balance (F1.as, F1.ab)
% are updated as well. Damping of the mass-balance feedback (CtrlVar.MassBalanceGeometryFeedbackDamping) is not
% implemented. The assembly needs these fields to be consistent with F1.h.
%
% On return RunInfo.Forward.uvhConverged is set. This is used by uvhRootFinding (which calls uvh2NotConvergent, and thus
% reduces the time step, if it is false), Ua2D and AdaptiveTimeStepping. The solution is considered to have converged if
% norm(f)^2 is less than Opt.ForceTolerance (the same criterion as the one used for rForce in SSTREAM_TransientImplicit).
% fsolve is stopped, by an output function, as soon as this is the case, and the exit flag returned by fsolve is not used.
%
% Thickness constraints:
%
% No active set is used inside this function, but it can be used from within uvhRootFinding (CtrlVar.ThicknessConstraints=true).
% The active nodes in BCs1.hPosNode are then treated as fixed values, as they are in the KKT system in SSTREAM_TransientImplicit,
% and on return the Lagrange multipliers l1.ubvb and l1.h are calculated from the converged residual, which is what
% ActiveSetUpdate needs. At the solution R+Aeq'*l=0 (the sign convention used in SSTREAM_TransientImplicit), and l is calculated as
% the minimum-norm least-squares solution of Aeq'*l=-R.
%
% The thickness penalty term in the assembly (CtrlVar.ThicknessPenalty=true) can be used with or without the active set.
%
% Only the SSTREAM flow approximation is supported.
%
% Options can be set through the (optional) structure CtrlVar.uvhFsolve. The fields, with defaults, are:
%
%   CtrlVar.uvhFsolve.Algorithm="trust-region-dogleg"  % or "trust-region" or "levenberg-marquardt"
%   CtrlVar.uvhFsolve.MaxIterations=50
%   CtrlVar.uvhFsolve.FunctionTolerance=[]             % [] : set to ForceTolerance^2. (fsolve declares success when norm(f)^2 < sqrt(FunctionTolerance))
%   CtrlVar.uvhFsolve.StepTolerance=1e-12
%   CtrlVar.uvhFsolve.OptimalityTolerance=1e-20        % fsolve's first-order optimality is absolute and depends on the scaling of J, so it is not used as the main stopping criterion
%   CtrlVar.uvhFsolve.ForceTolerance=[]                % [] : set to CtrlVar.uvhDesiredWorkAndForceTolerances(2). Convergence if norm(f)^2 is below this value
%   CtrlVar.uvhFsolve.Display="iter"                   % "off", "iter", "final"
%   CtrlVar.uvhFsolve.CheckGradients=false             % finite-difference check of J at the starting point (uses checkGradients). Very expensive: one residual evaluation per free variable. Only for tiny meshes.
%   CtrlVar.uvhFsolve.Scale=[]                         % see below
%
% Scaling of f:
%
%   Scale=[]         f is scaled with the norm of the projected external force vector, N'*Fext0, where Fext0 is the assembled
%                    force vector with ZeroFields=true (this is the same normalisation as used in SSTREAM_TransientImplicit).
%   Scale="initial"  f is scaled so that norm(f)=1 at the starting point.
%   Scale=number     f is multiplied with this number.
%
% The tolerances are therefore relative to the chosen scale.
%
% See also: PairwiseConstraintsNullSpace, uvh2D, uvhRootFinding, SSTREAM_TransientImplicit
%%

narginchk(8,8)

persistent WarnedAboutDamping

n=MUA.Nnodes;

if ~ismember(lower(string(CtrlVar.FlowApproximation)),["sstream","sstream-rho"])
    error("uvh_GaussNewton_fsolve:FlowApproximation","uvh_GaussNewton_fsolve is only implemented for the SSTREAM flow approximation.")
end

if CtrlVar.MassBalanceGeometryFeedback>0 && isfield(CtrlVar,"MassBalanceGeometryFeedbackDamping") && CtrlVar.MassBalanceGeometryFeedbackDamping~=0 && isempty(WarnedAboutDamping)
    warning("uvh_GaussNewton_fsolve:Damping","CtrlVar.MassBalanceGeometryFeedbackDamping is not implemented in uvh_GaussNewton_fsolve and is ignored.")
    WarnedAboutDamping=true;
end

%% Options

Opt.Algorithm="trust-region-dogleg";
Opt.MaxIterations=50;
Opt.FunctionTolerance=[];
Opt.StepTolerance=1e-12;
Opt.OptimalityTolerance=1e-20;
Opt.ForceTolerance=[];
Opt.Display="iter";
Opt.CheckGradients=false;
Opt.Scale=[];

if isfield(CtrlVar,"uvhFsolve") && isstruct(CtrlVar.uvhFsolve)
    fn=fieldnames(CtrlVar.uvhFsolve);
    for i=1:numel(fn)
        if isfield(Opt,fn{i})
            Opt.(fn{i})=CtrlVar.uvhFsolve.(fn{i});
        else
            warning("uvh_GaussNewton_fsolve:UnknownOption","Unknown option CtrlVar.uvhFsolve.%s ignored.",fn{i})
        end
    end
end

if isempty(Opt.ForceTolerance)
    if isfield(CtrlVar,"uvhDesiredWorkAndForceTolerances")
        Opt.ForceTolerance=CtrlVar.uvhDesiredWorkAndForceTolerances(2);   % same criterion as in SSTREAM_TransientImplicit
    else
        Opt.ForceTolerance=1e-15;
    end
end
% fsolve declares success if the sum of squares, norm(f)^2, is less than sqrt(FunctionTolerance)
if isempty(Opt.FunctionTolerance)
    Opt.FunctionTolerance=Opt.ForceTolerance^2;
end

%% Constraints

MLC=BCs2MLC(CtrlVar,MUA,BCs1);
if numel(l1.ubvb)~=numel(MLC.ubvbRhs) ; l1.ubvb=zeros(numel(MLC.ubvbRhs),1) ; end
if numel(l1.h)~=numel(MLC.hRhs) ; l1.h=zeros(numel(MLC.hRhs),1) ; end

% Here the constraints are taken directly from the MLC structure and not from AssembleLuvhSSTREAM. The latter can apply a
% mass-matrix scaling (if CtrlVar.LinFEbasis is true), and its output is meant for the KKT system and not for elimination.
Luv=MLC.ubvbL ; Lh=MLC.hL ;
if isempty(Luv) ; Luv=sparse(0,2*n) ; end
if isempty(Lh)  ; Lh=sparse(0,n) ; end
nu=size(Luv,1) ; nh=size(Lh,1) ;

Aeq=[Luv sparse(nu,n) ; sparse(nh,2*n) Lh];
beq=[MLC.ubvbRhs(:) ; MLC.hRhs(:)];

[N,xp,NInfo]=PairwiseConstraintsNullSpace(Aeq,beq);
ny=size(N,2);

%% Starting point

% As in SSTREAM_TransientImplicit: make sure that the geometrical fields are consistent with F1.h
[F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);

% The starting point is projected onto the set of feasible points.
x0=[F1.ub;F1.vb;F1.h];
y0=N.'*(x0-xp);

if CtrlVar.InfoLevelNonLinIt>=1
    fprintf("uvh_GaussNewton_fsolve: %i dofs, %i fixed, %i free (after elimination of %i fixed-value and %i tie constraints). |x0-x0feasible|=%g \n",...
        numel(x0),NInfo.nFixedDofs,ny,NInfo.nFixedRows,NInfo.nTieRows,norm(x0-(xp+N*y0)))
end

%% Solve

if ny>0

    f0=reducedFunction(y0,1,N,xp,UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

    if any(~isfinite(f0))
        error("uvh_GaussNewton_fsolve:NonFiniteResidual","The residual is not finite at the starting point.")
    end

    % scale
    normf0=norm(f0);
    scale=nan;

    if isnumeric(Opt.Scale) && ~isempty(Opt.Scale)

        scale=Opt.Scale;

    elseif isempty(Opt.Scale)

        CtrlVarZero=CtrlVar;
        CtrlVarZero.uvhMatrixAssembly.ZeroFields=true;
        CtrlVarZero.uvhMatrixAssembly.Ronly=true;
        [~,~,Fext0]=uvhAssembly(UserVar,RunInfo,CtrlVarZero,MUA,F0,F1,l1,BCs1);
        normFext=norm(N.'*Fext0);
        if normFext>0 && isfinite(normFext)
            scale=1/normFext;
        end

    elseif string(Opt.Scale)~="initial"

        error("uvh_GaussNewton_fsolve:Scale","CtrlVar.uvhFsolve.Scale must be empty, ""initial"", or a number.")

    end

    if isnan(scale)   % "initial", or fall back to this if external forces are zero
        if normf0>0
            scale=1/normf0;
        else
            scale=1;
        end
    end

    fun=@(y) reducedFunction(y,scale,N,xp,UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

    options=optimoptions("fsolve",...
        "Algorithm",Opt.Algorithm,...
        "SpecifyObjectiveGradient",true,...
        "Display",Opt.Display,...
        "MaxIterations",Opt.MaxIterations,...
        "FunctionTolerance",Opt.FunctionTolerance,...
        "StepTolerance",Opt.StepTolerance,...
        "OptimalityTolerance",Opt.OptimalityTolerance,...
        "OutputFcn",@(x,optimValues,state) StopWhenSolved(optimValues,state,Opt.ForceTolerance));

    if Opt.CheckGradients
        % fsolve itself no longer has a CheckGradients option. Use the stand-alone function checkGradients instead (MATLAB R2023b or later).
        [validJ,errJ]=checkGradients(fun,y0,Display="on"); %#ok<ASGLU>
        if ~validJ
            warning("uvh_GaussNewton_fsolve:Jacobian","The Jacobian does not agree with its finite-difference approximation.")
        end
    end

    tSolve=tic;
    [y,fval,exitflag,output]=fsolve(fun,y0,options);
    tSolve=toc(tSolve);

else

    % all dofs are fixed
    y=zeros(0,1); fval=zeros(0,1); exitflag=1; scale=1; normf0=0; tSolve=0;
    output.iterations=0; output.funcCount=0; output.message="All variables are fixed.";

end

%% Map back to the full set of variables

x=xp+N*y;

F1.ub=x(1:n);
F1.vb=x(n+1:2*n);
F1.h=x(2*n+1:3*n);

% As in SSTREAM_TransientImplicit
[F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow,F1.GF);
if CtrlVar.MassBalanceGeometryFeedback>0
    F1=UpdateMassBalance(UserVar,CtrlVar,MUA,F1);
end

%% Lagrange multipliers

% At the solution R+Aeq'*l=0, where R is the (unprojected) residual. This is the sign convention used in
% SSTREAM_TransientImplicit. The multipliers are needed by ActiveSetUpdate. The rows of Aeq are not necessarily linearly
% independent, hence lsqminnorm.

if size(Aeq,1)>0

    R=reducedFunction(x,1,speye(3*n),zeros(3*n,1),UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);
    lMult=-lsqminnorm(Aeq.',R);
    l1.ubvb=lMult(1:nu);
    l1.h=lMult(nu+1:end);

    if CtrlVar.InfoLevelNonLinIt>=10
        fprintf("uvh_GaussNewton_fsolve: |R+Aeq'*l|/|R|=%g \n",norm(R+Aeq.'*lMult)/norm(R))
    end

end

%% Convergence and RunInfo

relRes=norm(fval);
converged = relRes^2<=Opt.ForceTolerance && all(isfinite(x)) ;

RunInfo.Forward.uvhConverged=converged;

iStep=CtrlVar.CurrentRunStepNumber;
if numel(RunInfo.Forward.uvhIterations) < iStep
    RunInfo.Forward.uvhIterations=[RunInfo.Forward.uvhIterations;RunInfo.Forward.uvhIterations+NaN];
    RunInfo.Forward.uvhResidual=[RunInfo.Forward.uvhResidual;RunInfo.Forward.uvhResidual+NaN];
    RunInfo.Forward.uvhBackTrackSteps=[RunInfo.Forward.uvhBackTrackSteps;RunInfo.Forward.uvhBackTrackSteps+NaN];
end
RunInfo.Forward.uvhIterations(iStep)=output.iterations;
RunInfo.Forward.uvhResidual(iStep)=relRes^2;
RunInfo.Forward.uvhBackTrackSteps(iStep)=NaN;

if CtrlVar.InfoLevelNonLinIt>=1
    fprintf("uvh_GaussNewton_fsolve: converged=%i, exitflag=%i, iterations=%i, function evaluations=%i, |f|^2=%g (|f0|^2=%g), scale=%g, CPU=%gs \n",...
        converged,exitflag,output.iterations,output.funcCount,relRes^2,(normf0*scale)^2,scale,tSolve)
end

if ~converged
    warning("uvh_GaussNewton_fsolve:NoConvergence","uvh_GaussNewton_fsolve did not converge: exitflag=%i, |f|^2=%g (tolerance %g). %s",exitflag,relRes^2,Opt.ForceTolerance,strtrim(string(output.message)))
end

end


function [f,J]=reducedFunction(y,scale,N,xp,UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)

%%
% Projected and scaled residual,  f=scale*N'*R,  and its Jacobian,  J=scale*N'*K*N
%
% Before the assembly, the geometry (and possibly the mass balance) are updated to be consistent with the new iterate, in the
% same way as done in the Newton loop in SSTREAM_TransientImplicit.
%%

n=MUA.Nnodes;

x=xp+N*y;
F1.ub=x(1:n);
F1.vb=x(n+1:2*n);
F1.h=x(2*n+1:3*n);

% update s and b
if ~CtrlVar.ResetThicknessInNonLinLoop
    CtrlVar.ResetThicknessToMinThickness=0;
end
[F1.b,F1.s]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);

% update mass balance
if CtrlVar.MassBalanceGeometryFeedback>0
    F1=UpdateMassBalance(UserVar,CtrlVar,MUA,F1);
end

CtrlVar.uvhMatrixAssembly.ZeroFields=0;

if nargout==1

    CtrlVar.uvhMatrixAssembly.Ronly=1;
    [~,~,R]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);
    f=full(scale*(N.'*R));

else

    CtrlVar.uvhMatrixAssembly.Ronly=0;
    [~,~,R,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);
    f=full(scale*(N.'*R));
    J=scale*(N.'*K*N);

end

end


function F1=UpdateMassBalance(UserVar,CtrlVar,MUA,F1)

% The mass balance is evaluated at the end of the time step (as in SSTREAM_TransientImplicit)
tOld=CtrlVar.time;
CtrlVar.time=tOld+CtrlVar.dt;
F1.time=CtrlVar.time;

[~,F1]=GetMassBalance(UserVar,CtrlVar,MUA,F1);

F1.time=tOld;

end

function stop=StopWhenSolved(optimValues,state,ForceTolerance)

%%
% Output function for fsolve. Stops the iteration as soon as norm(f)^2 is less than ForceTolerance, which is the same type of
% criterion as used for rForce in SSTREAM_TransientImplicit.
%%

stop=false;

if state=="iter"
    stop = norm(optimValues.fval)^2 <= ForceTolerance ;
end

end
