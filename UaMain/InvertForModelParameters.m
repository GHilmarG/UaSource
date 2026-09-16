
function [UserVar,F,l,InvFinalValues,RunInfo]=...
    InvertForModelParameters(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

narginchk(11,11)

%%
%
%
% As shown in detail in the Ua compendium, an inversion for the parameters p, can be formulated as a constrained minimization
% problem:
%
% $$J(p)=I(q(p),p) + R(p) $$
%
% subject to
%
% $$F(q(p),p)=0$$
%
% To use a gradient-based optimization method, we need to be able to calculate the derivatives of J with respect to p.
%
% Using the Lagrange method, we form the extended functional
%
% $$\mathcal{L}(p)=I(q(p),p) + R(p) + \langle F(q(p),p) | \lambda \rangle $$
%
%
% We can then find the directional derivative of $J(p)$ with respect to $p$ in the direction $\phi$ as follows:
%
% 1) For some $p$, solve for $q$ using the forward model:
%
% $$F(q(p),p)=0$$
%
% 2) Then for this $q$, solve the linear adjoint problem:
%
% $$ \langle (d_q F)^* \lambda | \phi \rangle = - \langle  d_q J | \phi \rangle  $$
%
% 3) And then evaluate the total derivative with respect to $p$ as
%
% $$ d_p J = \langle (d_q F)^* | \lambda  \rangle + \partial_p J $$
%
% As an example, consider a $B$ inversion using momentum and mass conservation (for grounded ice where $h=s-B$):
%
%
% $$J(B)= I(v(B),B) + R(B) $$
%
% subject to
%
% $$P(v(B),B) = 0 $$   (momentum)
%
% $$M(B) = 0 $$   (mass)
%
% we could introduce Lagrange multipliers for both of these equations, but we can also use the fact that $\dot{h}$ can so
% easily be calculated from the mass conservation equation, and put the mass-conservation directly into the cost function $J$
%
% $$J(B)= \| v_c - v_m \| + \| \dot{h}_c - \dot{h}_m \|  + \|B_c - B_m \| + \langle P(v(B),B) | \lambda \rangle $$
%
% we can write this as
%
% $$J(B)= I_v + I_{\dot{h}}  + R(B) + \langle P(v(B),B) | \lambda \rangle $$
%
% were $I_v$ and $I_{\dot{h}}$ are misfit terms,  and $R$ a regularization term, and where we simply calculate evaluate
%
% $$\dot{h}_c=a - \nabla (v \, (s-B) ) $$
%
% and insert into the misfit term, i.e.
%
% $$I_{\dot{h}} = \| (a - \nabla (v \, (s-B) ) - \dot{h}_m \| $$
%
% The misfit term $I_{\dot{h}}$ is an explicit function of $v$ (i.e. the q variable).
%
% The (linear) adjoint problem now reads
%
% $$ \langle (d_v F)^* \lambda | \phi \rangle = - \langle  d_v I_v | \phi \rangle  - \langle  d_v I_{\dot{h}} | \phi \rangle
% $$
%
% The right-hand term is the derivative of $J=I+R$ with respect to $v$, but as $R$ does not depend on $v$, we only have those
% two terms involving $I$.
%
% The directional derivative of $J$ with respect to $B$ is calculated from
%
%
% $$ d_B J = \langle (d_B F)^* | \lambda \rangle + \partial_B J $$
%
% * $$B$$ inversion
%
%
% $$J(B) = \frac{1}{2}\left(d_{obs}-d_{modelled} \right)^{T}C_d^{-1}\left(d_{obs}-d_{modelled}\right) + \frac{1}{2}\left(B_{obs}-G_{obs}B\right)^{T}C_{B_{obs}}^{-1}\left(B_{obs}-G_{obs}B\right) + \frac{1}{2}\left(B-B_{prior}\right)^{T}Q\;\left(B-B_{prior}\right)$$
%
% Gauss-Newton system
%
% $$\left(J_f^{T}C_d^{-1}J_f + G_{obs}^{T}C_{B_{obs}}^{-1}G_{obs} + Q\;\right) \;\Delta B = J_f^{T}C_d^{-1}\left(d_{obs}-d_{modelled}\right) + G_{obs}^{T}C_{B_{obs}}^{-1}\left(B_{obs}-G_{obs}B\right) + Q\;\left(B_{prior}-B\right)$$
%
% $$J_f$$ is the directional derivative of the forward model
%
% $$d_{modelled} = f(B) $$
%
% that is
%
%
% $$D_{\delta B}f(B) = J_f\,\delta B = \lim_{\epsilon\to0}\frac{f(B+\epsilon\,\delta B)-f(B)}{\epsilon} $$
%
% $$\langle J_f\,\delta B,\, w\rangle_Y = \langle \delta B,\, J_f^{*}w\rangle_X$$
%
% where $X$ and $Y$ are infinite-dimensional Hilbert spaces. The adjoint PDE is derived from the (continuous) forward model
% before it is discretized.
%
%%



%%

if ~isfield(MUA,'M') || isempty(MUA.M)
    MUA.M=MassMatrix2D1dof(MUA);
end

if ~isfield(MUA,'Dxx') || isempty(MUA.Dxx)
    [MUA.Dxx,MUA.Dyy]=StiffnessMatrix2D1dof(MUA);
end
%%  Are we using the Riez-mapped gradient?
%
%
% The answer to that is: yes when using a gradient-based approach, ie where the Hessian is not build directly.
%
% But how the Riez-mapped gradient is introduced differs deepening on if the Ua in-build conjugate-gradient optimiser or the
% MATLAB fmincon optimizer is used:
%
% If:
%
%
%   CtrlVar.Inverse.MinimisationMethod="-UaOptimization-GradientBased-"
%
% then the metric matrix, G, does the mapping as
%
% $$ \nabla_{H^1} J =G^{-1} \nabla_{l^2} J $$
%
% For
%
%   CtrlVar.Inverse.MinimisationMethod="-MatlabOptimization-GradientBased-"
%
% the Riez-based gradient is introduced through a change of variables and the minimisation is on
%
% $$J(u)$$
%
% instead of
%
% $$J(p)$$
%
% where
%
% $$ u = R p $$
%
% and
%
% $$\nabla J = R^{-T} \nabla_{l^2} J $$
%
%
% where $R$ is the Cholesky factorisation of $G$ with
%
% $$G=R^{T} R $$
%
% This ensures that
%
%
%   fmincon
%
% with the lBFGS update, uses a consistent pair of $J$ and $\nabla J$. This is essential for the MATLAB optimizer to use the
% right slope $\nabla J \cdot d$ where $d$ is the search direction, and correct evaluation of any Armijo and curvature tests
% it may do internally, and for descent to be guaranteed. The reason why a non-descent direction might be selected by MATLAB
% fmincon without a warning if feeding it with
%
% $$\nabla J = R^{-T} \nabla_{l^2} J $$
%
% is subtle. fmincon with lBFGS update, selects a search direction as
%
% $$d=-M g_{L^2}$$
%
% where $M$ is its lBFGS approximation to the inverse Hessian. If we feed it with the Riez-mapped
% gradient we get the search direction
%
%
% $$d= -M G^{-1}  g_{l^2} $$
%
% but while both $M$ and $G$ are SPD, the produce $M G^{1}$ need not be! Therefore $d$ is not guaranteed to be a descent
% direction.
%
%
%%



if CtrlVar.Inverse.RieszMapGradient
    fprintf(" The optimisation will use a Riesz-mapped gradient.\n")
else
    fprintf(" The optimisation will not use a Riesz-mapped gradient.\n")
end

if CtrlVar.Inverse.CholeskyMappingOfCostFunctionAndGradient
    fprintf(" The optimisation will use a Cholesky mapping of cost function and gradient.\n")
else
    fprintf(" The optimisation will not use a Cholesky mapping of cost function and gradient.\n")
end


%% What inversions are being performed?
%  And make sure the Matern parameters are all correctly defined
[CtrlVar] = TikhonovToMaternMapping(CtrlVar,MUA);

%%
[isA,isB,isC] = isABC(CtrlVar);

if isB
    if ~isempty(Meas.Bobs)
        if isempty(Meas.BO)
            % I'm expecting this field to have been populated already (this should have been done in GetInputsForInverseRun.m) but if
            % is has not, do it here. 
            % Create the NBode2DataMap matrix O, and keep information about meas inside/outside mesh
            [Meas.BO,Meas.BInside,Meas.BEleID]=BuildNode2DataMap(CtrlVar,MUA,Meas.Bx,Meas.By) ;
        end
    end
end


%% Now build the metric/precision matrix once and keep in MUA. This is OK because in an inversion MUA never changes
%
% The metric matrix is here built from the blocks of the precision matrices, it plays a number of different roles. It is the
% precision matrix in the inversion, the metric matrix when using the trust-region approach, and the inner-product matrix
% when defining the inner product and the steepest descent direction.
%
% I'm storing here individual blocks as well as the block matrix (MetricMatrix). This is a bit of a waste of memory, but
% these are all sparse matrices. But this could be revisited at a later stage.
[MUA.G,MUA.QA,MUA.QB,MUA.QC]=BuildMetricMatrix(CtrlVar,MUA,Meas);
MUA.dG=decompositionUa(MUA.G);
[MUA.RG,flag,MUA.PRG]=chol(MUA.G) ;


%% Define inverse parameters and anonymous function returning objective function, directional derivative, and Hessian
%

% F should always be populated with the fields used at each stage of the inversion, and all the calculations are done using
% F. Start by populating F with the starting values
F=InvStartValues2F(CtrlVar,MUA,F,InvStartValues,Priors,Meas) ;
[F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);


F.GF=IceSheetIceShelves(CtrlVar,MUA,F.GF) ;

% p is the vector of the control variables, currently p=[A,B,C]
% with A, B or C here only being nonempty when inverted for,
% This mapping between A, B and C into the control variable is done by F2p

% Make sure initial point is feasible
F.AGlen=kk_proj(F.AGlen,F.AGlenmax,F.AGlenmin) ;
F.B=kk_proj(F.B,F.Bmax,F.Bmin) ;
F.C=kk_proj(F.C,F.Cmax,F.Cmin) ;

% Get an initial correct forward solve.
[UserVar,RunInfo,F,l]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l) ;

%%  Do the mapping from F to p :
% This does the mapping from the (control) variables that we are inverting for (one or more of logA, B, and logC), to the
% control parameter vector p. The parameter that we are inverting for are therefore all contained in the variable p. p0 is
% the starting value.
[p0,plb,pub]=F2p(CtrlVar,MUA,F);


CtrlVar.Inverse.ResetPersistentVariables=1;


% JGH: Returns the cost function (J), the gradient of the cost function with respect to p (dJdp), and the Hessian (ddJddp).
% The Hessian of the regularization term (R) can usually be calculated exactly, while the Hessian of the misfit/likelihood term
% (I), can not. However, one can come up with an educated guess for the Hessian of I with respect to C.

CtrlVar.JGH.CalcHessian=false;
[J0,dJdp,~,JGHouts,F]=JGH(p0,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint);

CtrlVar.Inverse.ResetPersistentVariables=0;
% The parameters passed in the anonymous function are those that exist at the time the anonymous function is created.



% Function handles are created to the functions calculating the cost function, J, the gradient, dJdp, and the Hessian. This
% is then passed to the optimization libraries.

% It is not easy to pass updated information to the cost function after it has been defined. The only updated variable in
% each call is p itself.  I need to decide at this stage if the Hessian will ever be needed.


if contains(CtrlVar.Inverse.MinimisationMethod,"Hessian")
    CtrlVar.JGH.CalcHessian=true;  % From now on I CAN use JGH to calculate the Hessian.
    % But JGH will only do so if the number of output arguments is ALSO 3 or greater.
    %
    % This means that the number of output arguments to JGH
    % (and therefore to func), controls if the Hessian is calculated or not.
else
    CtrlVar.JGH.CalcHessian=false; % From now on I can NOT use JGH to calculate the Hessian.
    %
    % This means that irrespective of the number of output arguments to JGH
    % (and therefore to func), the Hessian will never be calculated by JGH (or func).
end

func=@(p) JGH(p,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint);   % returns the cost (J), gradient (G) and Hessian (H)
% The Hessian output is used with the UaOptimisation toolbox, and when using the trust-region-reflective algorithm


% Somewhat annoyingly when using the interior-point algorithm, the MATLAB optimization toolbox wants the Hessian returned in
% a separate function, so I can't use JGH (!?). The function HessianABC is just a wrapper around JGH and returns the same
% Hessian as JGH.
%
% But when using the trust-region-reflective algorithm, the Hessian is returned as the third output to JGH and the Hfunc is
% not needed.
Hfunc=@(p,lambda) HessianABC(p,lambda,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint); % returns the Hessian (H) for the interior-point method


%%


Aineq=[];
bineq=[];


fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \n \n',J0,JGHouts.MisfitOuts.I,JGHouts.RegOuts.R,norm(dJdp))

dJdpTest=[];

%%

if CtrlVar.Inverse.TestAdjoint.isTrue
    %% The correctness of the gradient calculation can be tested by comparing it with a brute-force finite differences calculations.

    [J,dJdp,dJdpTest] = TestCorrectnessOfAdjointGradient(func,p0,MUA,CtrlVar,plb,pub);


else


    %%

    if contains(CtrlVar.Inverse.MinimisationMethod,"Ua")

        [p,UserVar,RunInfo]=UaOptimisation(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub);

    elseif contains(CtrlVar.Inverse.MinimisationMethod,"Matlab")

        clear fminconOutputFunction fminconHessianFcn fminuncOutfun
        [p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub,Hfunc,Aineq,bineq);

    else
        error("CaseNotFound")
    end

    % Here the final values from inversion, which are in the vector p, are copied across to the corresponding fields of F
    F=p2F(CtrlVar,MUA,p,F,Meas,Priors);

    % And a final additional call is made to get the cost function, J, and the gradient, dJdp, at the end of the optimization. In
    % principle, I guess it should be possible to get this information from the (external) optimization subroutine, but I don't
    % know how...
    CtrlVar.JGH.CalcHessian=false;  % Make sure that I don't calculate the Hessian again here as well.
    [J,dJdp,~,JGHouts,F]=JGH(p,plb,pub,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint);
    fprintf('\n +++++++++++ At end of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \n \n',J,JGHouts.MisfitOuts.I,JGHouts.RegOuts.R,norm(dJdp))


end

% Put RAa, RAs, RCa, RCs in InvFinalValues
InvFinalValues=Vars2InvValues(CtrlVar,F,InvStartValues,J,dJdp,JGHouts,RunInfo,dJdpTest);


end



%%%%%%%%%%%%%%%%%%  LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function [CtrlVar] = TikhonovToMaternMapping(CtrlVar,MUA)

[isA,isB,isC] = isABC(CtrlVar);

% make sure the Matern parameters are correct and can be used, even if the old Tikhonov approach is still being used
if CtrlVar.Inverse.Methodology=="-Tikhonov-"

    if isA
        [CtrlVar.Inverse.Matern.logAGlen.alpha,CtrlVar.Inverse.Matern.logAGlen.tau,CtrlVar.Inverse.Matern.logAGlen.kappa]=Tikhonov2MaternParameters(CtrlVar.Inverse.Regularize.logAGlen.ga,CtrlVar.Inverse.Regularize.logAGlen.gs,MUA.Area);
    end

    if isB
        [CtrlVar.Inverse.Matern.B.alpha,CtrlVar.Inverse.Matern.B.tau,CtrlVar.Inverse.Matern.B.kappa]=Tikhonov2MaternParameters(CtrlVar.Inverse.Regularize.B.ga,CtrlVar.Inverse.Regularize.B.gs,MUA.Area);
    end

    if isC
        [CtrlVar.Inverse.Matern.logC.alpha,CtrlVar.Inverse.Matern.logC.tau,CtrlVar.Inverse.Matern.logC.kappa]=Tikhonov2MaternParameters(CtrlVar.Inverse.Regularize.logC.ga,CtrlVar.Inverse.Regularize.logC.gs,MUA.Area);
    end


end

% Although the user defines alpha, kappa and tau, I might find it useful to work in terms of rho, sigma and nu.
if isA
    [CtrlVar.Inverse.Matern.logAGlen.rho,CtrlVar.Inverse.Matern.logAGlen.sigma,CtrlVar.Inverse.Matern.logAGlen.nu]=Matern_alpha_kappa_tau(CtrlVar.Inverse.Matern.logAGlen.alpha,CtrlVar.Inverse.Matern.logAGlen.kappa,CtrlVar.Inverse.Matern.logAGlen.tau);
end
if isB
    [CtrlVar.Inverse.Matern.B.rho,CtrlVar.Inverse.Matern.B.sigma,CtrlVar.Inverse.Matern.B.nu]=Matern_alpha_kappa_tau(CtrlVar.Inverse.Matern.B.alpha,CtrlVar.Inverse.Matern.B.kappa,CtrlVar.Inverse.Matern.B.tau);
end
if isC
    [CtrlVar.Inverse.Matern.logC.rho,CtrlVar.Inverse.Matern.logC.sigma,CtrlVar.Inverse.Matern.logC.nu]=Matern_alpha_kappa_tau(CtrlVar.Inverse.Matern.logC.alpha,CtrlVar.Inverse.Matern.logC.kappa,CtrlVar.Inverse.Matern.logC.tau);
end
end





%%
function [J,dJdp,dJdpTest] = TestCorrectnessOfAdjointGradient(func,p0,MUA,CtrlVar,plb,pub)
% Get the gradient using the adjoint method

[isA,isB,isC] = isABC(CtrlVar); 

[J,dJdp]=func(p0);

NA=MUA.Nnodes;

% Find the subset (iRange) in p, for which the brute-force gradient is to be calculated
if isempty(CtrlVar.Inverse.TestAdjoint.iRange)
    nTests=min(20,numel(p0));    % just test for some random nodes
    iRange=randi(MUA.Nnodes,nTests,1);
else
    iRange=CtrlVar.Inverse.TestAdjoint.iRange;
end

I=(iRange>=1) & (iRange <= MUA.Nnodes);  % Just in case the use sets some CtrlVar.Inverse.TestAdjoint.iRange outside the nodal values in Mesh
iRange=iRange(I);

% if the inversion is done for more than one field, then expand iRange accordingly.
nBlocks=isA+isB+isC;

switch nBlocks
    case 2
        iRange=[iRange(:);iRange(:)+NA];
    case 3
        iRange=[iRange(:);iRange(:)+NA;iRange(:)+2*NA];
end

% Gradient calculated using a brute-force finite difference approach
dJdpTest = CalcBruteForceGradient(func,p0,plb,pub,CtrlVar,iRange);

Diff=norm(dJdp(iRange)-dJdpTest(iRange))/norm(dJdp(iRange));
fprintf("Test Adjoint gradients: Normalized differences between adjoint gradient and FD: %g \n ",Diff)

fig_dJdpTest=FindOrCreateFigure("Test dJdp") ; clf(fig_dJdpTest)
plot(dJdp(iRange),dJdpTest(iRange),"or") ; axis equal ;
hold on ;
plot([min(dJdp(iRange)) max(dJdp(iRange))],[min(dJdp(iRange)) max(dJdp(iRange))],"--k")
ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off
xlabel("Adjoint $\partial J/\partial p$ ",Interpreter="latex")  ;
ylabel("Finite difference $\partial J/\partial p$",Interpreter="latex")
title("$\partial J/\partial p$",Interpreter="latex")
subtitle(sprintf("Normalized diff %g",Diff),Interpreter="latex")


drawnow
end
%%