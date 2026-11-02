

%%
% 
% *Release Notes* _October 2026_
%
% * Support for compiled uvh, uv and h  matrix assemblies. 
% 
% Users can now use a complied version of the uvh, uv and h assembly routines. 
%
% Details are provided in : |UaCompiledUserGuide.pdf| 
% 
% MEX files for 64-bit Windows with MATLAB2026b are included in the distribution. 
% 
% To use this option set:
%
%   CtrlVar.uvhMatrixAssembly.CompiledKernel=true; 
%   CtrlVar.uvMatrixAssembly.CompiledKernel=true; 
%   CtrlVar.hMatrixAssembly.CompiledKernel=true; 
%
% and also set:
%
%   CtrlVar.Parallel.uvAssembly.spmd.isOn=false;    % the compiled uv-assembly is parallel already
%   CtrlVar.Parallel.uvhAssembly.spmd.isOn=false;   % the compiled uvh-assembly is parallel already
%
% For other platforms than Windows, consult  |UaCompiledUserGuide.pdf| 
%
%  * Support for the PARDISO direct sparse solver:
% 
%  Users can now use the PARDISO sparse linear solver as an alternative to the inbuilt MATLAB linear solvers (UMFPACK for the
%  uvh and CHOLMOD for the uv system).
%
% To use this option set:
%
%   CtrlVar.LinearSolver.Type="PARDISO";  
% 
% The advantage of PARIDISO over the default MATLAB solvers is that PARDISO allows for reuse of the analysis phase of the
% solution. This is particularly interesting when a linear system with the same sparsity structure is solved repeatedly. While
% PARDISO is not faster than the MATLAB solvers, the MATLAB solvers do not allow for this reuse of the analysis phase. In
% practice this analysis phase is quite expensive and a big part of the overall solve. Using PARIDSO can therefore
% significantly speed up the solve (factor of 5 to 10) and lead to reduction in overall run times of maybe 30 to 50%,
% although the exact numbers are situation dependent. Information about the sparsity pattern is now stored as part of the MUA
% structure and recalculated when the mesh changes. Some re-analyses phases may be required when using dynamical boundary conditions, but
% this has been optimized within Úa and is now only occasionally required. 
% 
% Installation guide is found in 
% 
%   UaPardisoGuide.pdf
%
%
%
% * Grounding-line plots now based on $h-h_{gl}$, ie thickness minus grounding line thickness, rather than the flotation mask.
% Thanks to Deborah Rhee for this suggestion. The grounding lines may now look a bit smoother. (Note that the grounding line
% position itself is never used in the code when finding the solutions.)
%
% * Default for the positive-thickness constraints: active set on, and thickness penalty off.
%
%   CtrlVar.ThicknessConstraints=true;
%   CtrlVar.ThicknessPenalty=false;
%
%
% * Active-set method:
% 
% Various changes were made to the active-set iteration to make it less prone to cycling and unnecessary and temporary
% release and activation of constraints. The active-set is not considerably more likely to converge in 1 to 3 iterations and
% often no updates are required. 
%
%
% * Explicit estimate:
%
% Explicit estimate is now, by default, only done for h. Boundary behavior related to level-set and min ice thickness are
% more robust. 
%
%
% * Reactions of the thickness constraints as an output:
%
% The equivalent local mass-balance rate of the thickness constraints is now provided as an output field, F.ah (units of distance per
% time, eg m/yr), positive where the constraints add ice and negative where they remove ice. It is calculated at the end of each time
% step from the discrete reactions of the constraints, and includes all thickness constraints: the active set (h>=ThickMin), the
% level-set constraints downstream of calving fronts, and user-defined thickness boundary conditions. In addition, the total rates at
% which mass is added and removed by the constraints (units of mass per time, eg kg/yr) are stored for each time step in
% RunInfo.Forward.ThicknessConstraintsMassAdded and RunInfo.Forward.ThicknessConstraintsMassRemoved. Downstream of a calving front, for
% example, the mass removed is the calving flux. The mass budget of each time step closes exactly, ie the change in ice mass is equal to
% the integrated mass balance, minus the net outflow across the boundary of the domain, plus the mass added and minus the mass removed by
% the constraints (tested to round-off in the Thule example). 
%
%
% * Local mesh refinement with CtrlVar.MeshSizeMax=NaN (2 November 2026):
%
% NewDesiredEleSizesAndElementsToRefineOrCoarsen2.m: with CtrlVar.MeshSizeMax=NaN (as in the Thule example), nodes outside the refinement
% ranges (CtrlVar.MeshAdapt.GLrange, CFrange, ...) kept a NaN target size. The threshold for local refinement is taken from the sorted size
% ratios at the fraction CtrlVar.LocalAdaptMeshRatio of all elements; since sort places NaN last, the threshold itself became NaN whenever
% fewer than that fraction of the elements had a finite ratio, and then no element was refined or coarsened, silently. (Thule: the 50 km
% mesh was refined, the 25 km mesh was not.) Three changes:
%
%   - The fallback that sets all target sizes to CtrlVar.MeshSize when no criterion applied (all target sizes equal to MeshSizeMax) now also
%     recognises a NaN MeshSizeMax (NaN==NaN is false).
%   - Nodes still without a finite target size keep their current size (size ratio 1: neither refined nor coarsened); the function no
%     longer returns NaN target sizes.
%   - The refinement and coarsening thresholds are taken from the finite ratios only. With all ratios finite they are the same as before.
%
% Thule: the 50 km mesh is refined exactly as before (12052 elements, 6081 nodes); the 25 km mesh is now refined with MeshSizeMax=NaN, and
% identically to a run with MeshSizeMax=MeshSize (29810 elements, 15005 nodes).
%
% * Inverse problem with the element-wise kH (1 November 2026):
%
% The gradient and Hessian calculations of the inverse problem now use the same grounding-line mask at the integration points as the
% forward uv assembly, ie the element-wise kH (CtrlVar.kHElementAlpha>0, default 0.02), so that they are consistent with it:
%
%   - IntegrationPointkH.m (new): kH for the grounding-line mask at the integration points, F.kHElement if set, otherwise calculated
%     from F (ElementWisekH.m), or the scalar CtrlVar.kH if CtrlVar.kHElementAlpha=0. Used by the uv assemblies (vectorised and compiled),
%     the compiled uvh assembly, and the routines below.
%   - Gradients and sensitivities: dIdCq.m, dIdBqGeneral.m, dFuvdC.m, dFuvdBGeneral.m. Hessian terms: FCC.m, FCuv.m, Fqq.m, FBC.m, FBB.m,
%     FBuv.m and WeertmanSecondOrderDerivatives.m (new optional input kHe). The water-depth mask (HeHint, as Hposint in the uv assembly)
%     and the nodal geometry (dGeometrydB.m, bFunc.m, Calc_bs_From_hBS.m) keep the scalar kH. dIdAq.m does not involve the mask.
%   - HeavisideApprox.m and DiracDelta.m: k.*(x-x0) instead of k*(x-x0), so that k can be an array of the size of x (identical for
%     scalar k).
%   - InvertForModelParameters.m: kH_e is evaluated once, from the thickness and bed at the start of the inversion, and kept fixed during
%     the inversion (it is not differentiable with respect to B; in B inversions it is not updated as B, and with it h over grounded areas,
%     changes). After the inversion it is recalculated from the final geometry before the next uv solve. Inversions previously used the
%     scalar kH (CtrlVar.kHElementAlpha was set to zero).
%   - calcVerticalSurfaceVelocity.m: element-wise kH, from its nodal input fields.
%   - TestInverseHessianEndToEnd.m: new options kHElementAlpha and kH; kH_e is kept fixed during the test, as during an inversion.
%
% Since HeavisideApprox.m and DiracDelta.m are compiled into the MEX kernels, all MEX variants were rebuilt.
%
% Tests (TestInverseHessianEndToEnd.m, synthetic problem, cases A, C, A+C with velocity data, B and A+B+C with velocity and dh/dt data):
% alpha=0: 5 of 5 passed. alpha=0.1 (kH_e<kH in 130 of 400 elements): 5 of 5 passed, gradients and Hessian with second-order convergence.
% alpha=0.2 (kH_e<kH in 367 of 400 elements): with the default steps the A and C cases fail, but with smaller steps the adjoint gradient
% agrees with finite differences to 1e-7 and the Hessian converges with second order: the cost function of this synthetic problem then
% varies strongly in curvature over changes of about 0.1% in A or C, which the default steps do not resolve. (With a scalar kH of
% comparable size, 0.26, all cases pass with the default steps; the origin of this non-linearity is not yet understood.)
%
% * Element-wise kH is the default, CtrlVar.kHElementAlpha=0.02 (1 November 2026):
%
% Cause of the non-convergence of the uvh Newton iteration in strongly transient phases with significant grounding-line migration (Thule
% example): the grounding-line mask GF=HeavisideApprox(kH,h-hf) is evaluated at the integration points. If the range of h-hf over a
% grounding-line element is much larger than the transition width 2/kH (Thule, 50 km mesh: about 1300 m against 2 m for kH=1), the
% integration points practically never lie within the transition, and the mask is effectively a step function of the thickness at the
% integration points. The residual is then nearly piecewise constant, with steep (smooth, but over about 1 m of thickness) transitions
% where an integration point passes through the transition. When a few integration points lie within the transition, the discrete
% problem may have no root near the iterate, only a local minimum of the residual norm (shown with fsolve: first-order optimality to
% 5e-6 at |f|^2=2e-7). Newton, kH continuation, linear rheology in ice-free areas, a regularised Jacobian and fsolve all stalled at the
% same residual level. The Jacobian is exact (Taylor test), the linear solves are exact, and there is no round-off floor.
%
% With the element-wise kH (ElementWisekH.m, see the entries below), the uvh solves in the Thule example converge without time-step
% reductions at all time steps up to dt=0.1 on all meshes tested (uniform 50, 25, 12.5 km; refined to 5.8, 2.9, 1.4 km at the grounding
% line), where the scalar kH=1 needs reductions on all of them; Newton iterations are 20-35% lower, and in a mesh-convergence test the
% element-wise and the scalar kH converge to the same solution (0.1-0.2 m rms difference on the finest meshes). In Greenland (real A and
% C, MAR mass balance, dt=0.2-1) both converge; the element-wise kH needs 4-22% fewer Newton iterations, with the same thickness accuracy
% and a systematic speed difference of about 1.75 m/yr rms, at a few grounding-line nodes.
%
% CtrlVar.kHElementAlpha=0.02 is therefore the default. This changes results slightly near grounding lines compared with earlier versions.
% Set CtrlVar.kHElementAlpha=0 to use the scalar kH everywhere, as before. A grounding-line fallback (the element-wise kH only in time
% steps in which the solve with the scalar kH stagnates, CtrlVar.uvhGLFallback), implemented earlier the same day, has been removed again:
% with the element-wise kH as the default it is not needed. If a uvh solve does not converge, the time step is reduced as before.
%
% * Compiled uv and uvh kernels with the element-wise kH (1 November 2026):
%
% The compiled uv and uvh kernels (uvElementKernel.m, uvhElementKernel.m and their entry points uvKernelFull/Ronly.m and
% uvhKernelFull/Ronly.m) take the element-wise kH as a new last input, kHv (Nele x 1), used for the grounding-line mask at the
% integration points (the nodal mask and Hposint keep the scalar kH). uvMatrixAssemblyCompiled.m and uvhMatrixAssemblyCompiled.m pass
% kHv as the vectorised assemblies do (from F.kHElement or F0.kHElement if set, otherwise from ElementWisekH.m; kH for all elements if
% CtrlVar.kHElementAlpha=0, which leaves the arithmetic unchanged). CompiledAssemblyCheck.m no longer excludes CtrlVar.kHElementAlpha>0.
% All eight uv and uvh MEX variants rebuilt (BuildUaCompiled.m). TestUaCompiled.m: new tests testUvhElementWisekH and testUvElementWisekH
% (compiled against vectorised assembly, kH_e calculated and taken from the field; the element-wise kH changes the residual).
% Regression tests: 90 of 90 passed.
%
% Run times with the element-wise kH (alpha=0.02), compiled against vectorised: Thule (refined mesh, 42706 nodes, dt=0.1) 55.5 s
% against 139.5 s (results identical to round-off), Greenland (56070 nodes, dt=1, t=0-6) 138.3 s against 476.7 s. Both are now faster
% than with the scalar kH (113.0 s and 169.4 s), owing to fewer Newton iterations and no time-step reductions.
%
% * Element-wise kH in the uv solver, computed once per time step (1 November 2026):
%
% The element-wise kH (CtrlVar.kHElementAlpha>0) is now defined in one place, ElementWisekH.m (from a state F) and ElementWisekHFromNodal.m
% (from the nodal values of h-hf of each element), and used by both the uvh and the uv assemblies:
%
%   - At the start of each time step (uvhRootFinding.m, uvhSemiImplicit.m), kH_e is evaluated from the state at the start of the time step
%     and stored in the new element field F.kHElement (Nele x 1; [] if CtrlVar.kHElementAlpha=0). The uv solves of -uv-h- use it, so that
%     -uvh- and -uv-h- solve the same discrete problem (Thule, alpha=0.02: the two schemes agree at t=1.5 to within 5 mm rms in thickness).
%   - uvMatrixAssemblySSTREAM.m: the grounding-line mask at the integration points uses kH_e (from F.kHElement if set, otherwise from the
%     current thickness, eg in diagnostic uv solves). The uv Jacobian remains consistent (Taylor test: second order).
%   - The nodal mask (Calc_bs_From_hBS.m) keeps the scalar kH.
%   - Inversions use the scalar kH (InvertForModelParameters.m sets CtrlVar.kHElementAlpha=0 with a message): the element-wise kH is not yet
%     implemented in the adjoint and gradient calculations.
%   - The compiled uv and uvh kernels do not yet support the element-wise kH (the vectorised assembly is used, CompiledAssemblyCheck.m), nor
%     does the spmd uv assembly (the sequential assembly is used, uvMatrixAssembly.m).
%   - F.kHElement is reset on mesh changes (MapFbetweenMeshes.m) and before diagnostic uv solves in Ua2D.m.
%   - CtrlVarValidityCheck.m: CtrlVar.kHElementAlpha must be finite and non-negative; a warning is given if it is larger than 0.5.
%
% With CtrlVar.kHElementAlpha=0 results are unchanged (regression tests, Greenland -uvh- and -uv-h- bit-identical).
%
% Tests of the element-wise kH used throughout (alpha=0.02), against the scalar kH=1:
%   - Thule (fixed calving front), mesh convergence on uniform (50, 25, 12.5 km) and refined meshes (5.8, 2.9, 1.4 km at the grounding line):
%     the two converge to the same solution; on the finest meshes they differ by 0.1-0.2 m rms in thickness, against discretisation errors
%     of 60-120 m rms on the coarser meshes. With dt=0.1 the scalar kH needs time-step reductions on all meshes, the element-wise kH none.
%     Newton iterations per time step are 20-35% lower with the element-wise kH.
%   - Greenland (real A and C, MAR mass balance, nx150), dt=0.2-1, t=0-6: all converge with both; 4-22% fewer Newton iterations with the
%     element-wise kH; thickness equally accurate (rms 0.4/0.1/0.03 m against dt=0.05 for dt=1/0.5/0.2); a systematic speed difference of
%     about 1.75 m/yr rms, concentrated at a few grounding-line nodes.
%
% * Convergence criteria of the uvh, uv, h and level-set solvers depend on rForce only, with the same logic (1 November 2026):
%
% uvhResidualsCriteria.m, used by the implicit uvh solver and for the outer uv-h iteration, now decides on the force residual rForce alone:
%
%   CtrlVar.uvhDesiredForceTolerance=1e-15;     % convergence
%   CtrlVar.uvhAcceptableForceTolerance=1e-6;   % acceptance of a stagnated solve, or of one ended by a short backtracking step
%
% Previously, the work residual rWork also entered: CtrlVar.uvhAcceptableWorkOrForceTolerances=[1 1e-8] required rWork<1 or rForce<1e-8
% for acceptance, in addition to rForce<1e-6. rWork contains the Newton step, which becomes very large at nearly singular modes of the
% Jacobian (eg weakly determined velocities of floating, nearly ice-free areas behind a calving front), and is returned as NaN (and then
% set to zero) by the line search whenever the accepted step is not a Newton step. Whether a stalled solve was accepted, or the time step
% halved, therefore depended on which line-search branch produced the last step (Thule example: a solve stalled at rForce=3e-7 after a
% Newton step was rejected and the time step halved; one stalled at rForce=5e-7 after a fallback step was accepted, and its solution was
% as accurate as that obtained with two half-steps). The acceptable tolerance 1e-6 is now the same as CtrlVar.uvhStalledForceTolerance in
% uvhRootFinding.m. The fields CtrlVar.uvhDesiredWorkAndForceTolerances, CtrlVar.uvhDesiredWorkOrForceTolerances,
% CtrlVar.uvhAcceptableWorkAndForceTolerances and CtrlVar.uvhAcceptableWorkOrForceTolerances are no longer used; if a run sets them, the
% force entries of the *WorkAndForceTolerances fields are taken over, with a warning (CtrlVarValidityCheck.m). The corresponding criteria
% for the h and level-set equations were changed in the same way, see below.
%
% The same change for the h equation (MassContinuityEquationNewtonRaphson.m, SSHEET_TransientImplicit.m), the level-set equation
% (LevelSetEquationNewtonRaphson.m) and the uv equations (uvResidualsCriteria.m, SSTREAM2dNR2.m): the criteria depend on rForce only, with
% new tolerances CtrlVar.hDesiredForceTolerance=1e-10, CtrlVar.hAcceptableForceTolerance=1e-6, CtrlVar.LSFDesiredForceTolerance=1e-12,
% CtrlVar.LSFAcceptableForceTolerance=1e-10, CtrlVar.uvDesiredForceTolerance=1e-15 and CtrlVar.uvAcceptableForceTolerance=1e-15, and the
% previous fields taken over with a warning if set. All solvers now follow the same logic as the uvh solver: the desired tolerance is
% required while the Newton iteration makes normal progress, and the acceptable tolerance is used if it makes poor progress (short or
% failed backtracking step, stagnation). In MassContinuityEquationNewtonRaphson.m and LevelSetEquationNewtonRaphson.m this was the other way
% round (the stricter tolerance with poor progress, contrary to the comments in the code). In the h solver, the hard-wired test for
% continuing while still converging quadratically uses rForce>1e-15 instead of rWork>1e-20. The uv solver (SSTREAM2dNR2.m) now, as the uvh
% solver, accepts a solution after failed backtracking if rForce is below the acceptable tolerance (previously always non-convergence),
% and detects stagnation (r/r0>0.98 in three successive iterations; previously the iteration continued to CtrlVar.NRitmax). The uv
% acceptable tolerance is kept at 1e-15 (unlike uvh), since uv solves are also used in inverse runs, where the adjoint requires a converged
% forward solution.
%
% Tests: h equation, replaying all 936 convergence checks of the h solves in three uv-h runs (Greenland with uniform and with inverted A
% and C, Thule), the new criteria make identical decisions; Greenland -uvh- and -uv-h- runs (dt=0.1) are bit-identical. Level-set equation:
% after the first Newton iteration rForce is about 1e-13; whether a second iteration followed was previously decided by rWork alone (in
% 177 of 206 solves of a Thule run with a moving calving front). Now all level-set solves converge after one iteration. Thule example,
% t=0 to 0.1: level-set iterations 77 -> 51; time steps, uvh Newton iterations and uvh residuals unchanged.
%
%
% * Active-set method: cumulative release criterion for thickness constraints (31 October 2026):
%
% A node constrained to the minimum thickness (h>=ThickMin) was released only if, without the constraint, the ice there would thicken
% by more than ReleaseAlpha*ThickMin within one time step. As a rate, this threshold is ReleaseAlpha*ThickMin/dt, so that slow
% accumulation on ice-free ground was suppressed, the more so the smaller the time step: with ReleaseAlpha=1e-3 and ThickMin=1 m, all
% accumulation below 1 m/yr was suppressed at dt=0.001. (Greenland, nx150, MAR mass balance: 84% of the ice-free nodes with positive
% mass balance were held at ThickMin at dt=0.001, against 54% at dt=0.01.)
%
% The default is now a cumulative criterion, CtrlVar.ActiveSet.ReleaseRule="cumulative": for each constrained node, the net thickening
% that the constraint has suppressed since the node became constrained is accumulated over time steps, as the new nodal field
% F.hSuppressedThickening (S, units of distance; S <- max(0, S-ah*dt) once per accepted time step, with ah the equivalent mass-balance
% rate of the constraint, see UpdateSuppressedThickening.m), and the node is released once max(0,S-ah*dt) > ReleaseAlpha*ThickMin. This
% does not depend on the time step: accumulation, however small, always leads to the release of the node eventually (after the ice
% would have grown by ReleaseAlpha*ThickMin, ie 1 mm with the defaults), while round-off contributions of alternating sign, the reason
% for ReleaseAlpha>0, do not accumulate. S is capped at the release threshold, is saved in restart files, and is mapped onto new meshes
% (MapFbetweenMeshes.m). The default value of ReleaseAlpha (1e-3) is unchanged; its exact value is much less important now. The previous
% rule can be selected with CtrlVar.ActiveSet.ReleaseRule="per-step".
%
% Greenland, nx150, MAR mass balance, -uvh-: with the cumulative rule, 42.5% of the ice-free nodes with positive mass balance are held at
% ThickMin at both dt=0.01 and dt=0.001 (the remainder lose more ice by outflow than they gain), and the ice thickness at t=0.5 differs
% by at most 7 mm between the two time steps, and by at most 4 mm from runs without any release threshold (ReleaseAlpha=0), which,
% however, release many nodes in ablation areas on round-off (25% of them, against 8% with the cumulative rule) at a much higher cost.
%
%
% * Active-set method: predictive activation and release threshold, new defaults (14 October 2026):
%
% Two changes to the active-set method are now the default:
%
%   CtrlVar.ActiveSet.PredictiveActivation=true;   % (previously not available)
%   CtrlVar.ActiveSet.ReleaseAlpha=1e-3;            % (previously 0)
%
% With predictive activation, nodes where the explicit estimate of the ice thickness at the end of the time step is at or below ThickMin
% (while the thickness at the start of the time step is above ThickMin) are added to the active set at the start of the time step.
% Previously, such nodes were only activated after the first uvh solve of the time step, which then always required a further solve. An
% optional margin can be set through CtrlVar.ActiveSet.PredictiveActivationMargin (default 0).
%
% With ReleaseAlpha>0, a constrained node is only released if the ice there would, without the constraint, thicken by more than
% ReleaseAlpha*ThickMin within the time step. With ReleaseAlpha=0, nodes were released as soon as their multiplier had the wrong sign. In
% some applications (eg Greenland) many of these releases are of nodes in ice-free areas whose multipliers are only negative at round-off
% level (relative thickness gains of order 1e-8), and many of these nodes immediately fell below ThickMin again, requiring further
% active-set iterations. 
%
%
% * The default value of CtrlVar.ATSTargetIterations has been changed from 4 to 7 (10 October 2026). Runs that do not set this parameter
% will therefore generally use larger time steps than before. To obtain the previous behaviour, set
% CtrlVar.ATSTargetIterations=4.
%
%
% * Estimate of the time-discretisation error of the ice thickness (10 October 2026):
%
% If CtrlVar.TimeDiscretisationErrorEstimate.Use=true, the local time-discretisation error of h is estimated in each time step
% (TimeDiscretisationErrorEstimate.m). This is a predictor-corrector estimate, comparing the second-order explicit estimate of h with the
% computed (implicit, theta=0.5) solution, and is only available for theta=0.5. The estimate is in units of h (ie m), and is the error
% committed in that time step. It is stored as a nodal field in F.hTimeDiscretisationErrorEstimate, and norms of it are stored as time series
% in RunInfo.Forward (hTimeDiscretisationErrorRMS, ..Max, ..Volume, ..ValidArea, etc.). Norms are calculated as sqrt(x'*M*x/A) with the
% consistent mass matrix M, where A is the area over which the estimate is defined.  Optionally, the estimate can be accumulated over a run
% (CtrlVar.TimeDiscretisationErrorEstimate.Accumulate=true, see F.hTimeDiscretisationErrorAccumulated). The accumulated estimate is an
% "error budget", ie the sum of the local errors, and not the global error, as errors are transported and may decay.
%
%
% AdaptiveTimeStepping.m is now a driver that selects between two approaches, depending on CtrlVar.AdaptiveTimeSteppingMethod:
%
%   "iteration-based" : the time step is based on the number of non-linear iterations (AdaptiveTimeSteppingIterationBased.m). This is the
%                       previous approach, unchanged, and is the default.
%   "error-estimate"  : in addition based on the estimated time-discretisation error of h (AdaptiveTimeSteppingWithErrorEstimate.m).
%
% With "error-estimate", the time step is chosen such that the estimated local error of h per unit time, scaled by atol+rtol*|h|, is about
% one, ie dt_new = Safety*dt*q^(-1/2), where q is the scaled error per unit time in the last time step (ATSTimeDiscretisationErrorLimit.m).
% The time step may then also be increased up to this accuracy-based time step, provided the convergence of the non-linear iterations is
% acceptable (number of iterations at most CtrlVar.ATSTargetIterations over the last few time steps); otherwise the minimum of the
% accuracy-based time step and the time step based on the non-linear iterations is used. Integrated over a run of length T, the accumulated
% local error is then approximately bounded by T*atol (for rtol=0), and the actual error at the end of a run decreases with atol
% (approximately as atol^0.7 in tests). The time-discretisation error estimate is switched on automatically. Settings, with default values:
%
%   CtrlVar.AdaptiveTimeSteppingMethod="iteration-based";      % or "error-estimate"
%   CtrlVar.ATSTimeDiscretisationError.atol=0.01;               % m/yr; values of about 0.0005 to 0.003 have given good results in tests
%   CtrlVar.ATSTimeDiscretisationError.rtol=0;                  % 1/yr
%   CtrlVar.ATSTimeDiscretisationError.StartTime=-inf;          % the accuracy-based time step is only used from this time onwards
%   CtrlVar.ATSTimeDiscretisationError.Safety=0.85; ...MaxIncrease=2; ...MaxDecrease=5; ...MinValidArea=0.1;
%   CtrlVar.ATSTimeDiscretisationError.GovernIncreases=true; ...ConvergenceCheckSteps=2;
%
% Recommendations: With "error-estimate", CtrlVar.ATSTargetIterations can be set to about 7 to 10, as the accuracy is then controlled
% separately.
%
%
% *Release Notes* _September 2026_
%
% * The conjugate gradient method used in  Ua optimisation has been completely rewritten, and is now competitive with the MATLAB lBFGS
% approach. This is the approach used when setting:
% 
%   CtrlVar.Inverse.MinimisationMethod="-UaOptimisation-GradientBased-" ;
% 
% Using this approach, the gradients are Sobolov gradients with respect to a metric defined by the explicit Hessians of the
% regularization part. Box constraints are implemented using two-sided variable transform.
% 
% Selecting inversion algorithms, e.g. gradient-based versus Hessian-based has been streamlined and
% simplified. This may possibly require some changes in input files. 
% Hessian solver now calculates the exact and full Hessian using the Direct-Adjoint method. This is very memory intensive and
% the Hessian option can only be used for small meshes (10,000 nodes or less). 
%
% UaHessian inverse-solver uses trust-region approach with a 2D exact sub-space minimisation.
%
%
% *  UK spelling used for `Optimisation' 
%
% *Release Notes* _August 2026_
%
%
% * For the -uvh- solve, two normalization of the residuals are now available, selected by
%
%    CtrlVar.uvhResidualNormalisation = "pooled"    (older approach)
%                                     = "blockwise" (new default approach)
%
% The new blockwise normalization option is selected by setting 
%
%    CtrlVar.uvhResidualNormalisation="blockwise";
%
% It scales the uv and the h blocks of the uvh residual vector separately. This should give a better measure of the residual and
% be less affected by relative differences in the uv and h residual blocks. 
%
% * There is a new "softplus" positive ice-thickness penalty formulation. It can be used in connection with the active-set iteration
% enforcing the min ice thickness constraint. It adds an implicit penalty term to the uvh system. The penalty term is a
% "soft" i.e. smooth function of the thickness violation. It has two parameters:
%
%   CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.K
%   CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.l
%
% $K$ is the slope of the linear term and l is a smoothness parameter. Note that $l$ has the units of ice thickness and should be
% selected to be somewhat smaller than CtrlVar.ThickMin
%
% Typical values might be
%
%   CtrlVar.ThicknessPenalty=true;
%   CtrlVar.ThicknessPenaltyMassBalanceFeedbackFunction="softplus";
%   CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.K=100;
%   CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.l=0.1*CtrlVar.ThickMin ;
%
% (Note that internally K is divided by the time step, dt. As a result added mass amplitude per uvh Newton iteration is independent of dt.)
% (No longer the case since 9 October 2026, K is again a rate. See the October 2026 release notes.)
%
% *Release Notes* _July 2026_
%
% An option added to use Matérn covariance matrices in an inversion. The resulting precision matrices are either:
%
%
% $$ Q_{\alpha=1} = \tau^2\left(\kappa^2 M + D\right) $$
%
% $$Q_{\alpha=2} = \tau^2 \; (\kappa^2 M + D)\; \tilde{M}^{-1} \; (\kappa^2 M + D)$$
%
% or
%
% $$Q_{\alpha=3} = \tau^2\left(\kappa^2M+D\right)\,\tilde{M}^{-1}\left(\kappa^2M+D\right)\,\tilde{M}^{-1}\left(\kappa^2M+D\right)$$
%
% where $M$ is the mass matrix, $\tilde{M}$ the lumped mass matrix, $D$ the stiffness matrix, and $\alpha$, $\kappa$ and $\tau$ are (hyper) parameters.  This option
% is activated by setting:
%
%  CtrlVar.Inverse.Methodology="-Matern-" ;
%
% and there is extensive description of this approach in the UaCompendium. The Matérn covariance function is:
%
% $$\left(C_M\right)_{ij} = C(r_{ij}) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}\left(\kappa r_{ij}\right)^{\nu}\, K_{\nu}\!\left(\kappa r_{ij}\right) $$
%
% where
%
% $$ r_{ij}=\|\mathbf{x}_i-\mathbf{y}_j\| $$
%
% Previously only the $\alpha=1$ option was implemented, through
%
%   CtrlVar.Inverse.Methodology="-Tikhonov-" 
%  
% which still is the default option. 
%
% When using the "-Tikhonov-" option, the parameters are entered as $\gamma_s$ and $\gamma_a$ in
%
% $$ Q_{\alpha=1} =  \gamma_a M + \gamma_s D  $$
%
% When using the new "-Matern-" option for $\alpha=1$ you can get the same behavior by selecting 
% $\tau$ and $\kappa$ in terms of $\gamma_s$ and $\gamma_a$ for $\alpha=1$. You can find the formulas in the UaCompendium or
% use the function: 
%
%   [alphaMatern,tauMatern,kappaMatern,sigma2Matern,nuMatern,rhoMatern]=Tikhonov2MaternParameters(ga,gs,Area)
% 
% Note that the relationship depends on the area of the computational domain. 
%
% *Release Notes* _June 2026_
%
% * An option added to test if prescribed boundary conditions (as prescribed by the user in DefineBoundaryConditions.m) are
% indeed linearly independent, and if not, uses a row-subselection to find the rows which are maximally independent. 
%
%
% Internally, the boundary conditions are represented as a constraint matrix 
%
% $$A_{\mathrm{eq}} \; x = b \quad ,  \quad A_{\mathrm{eq}} \in R^{m\times n} $$
%
% where $m$ is the number of constraints, and $n$ the number of degrees of freedom.
% 
% The matrix $A_{\mathrm{eq}}$ should have full row rank. It is really up to the user to make sure the BCs are in this sense consistent.
% However, it can be a bit tricky to ensure that this is the case for complicated BCs involving multiple links/ties etc.
% 
% By setting
%
%   CtrlVar.BCsRowSubsetSelection=true; 
%
% the matrix $A_{\mathrm{eq}}$ will be modified using a row-subset selection algorithm. This does involve a QR decomposition of $A_{\mathrm{eq}}$, but since $A_{\mathrm{eq}}$ is
% usually quite small, (few rows) this is fast. However, the best approach is for the user to make sure that this is not
% required in the first place by ensuring that the boundary conditions contain no redundancies. 
%
%
%
%
% * The default KKT solver for symmetrical matrices is now the Null-Space method. Previously the default was   
% 
%   CtrlVar.SymmSolver="EliminateBCsSolveSystemDirectly";
% 
% and this was used if the user used the default option
%
%    CtrlVar.SymmSolver='Auto'; 
% 
% This symmetrical solver preserves the symmetry of the reduced KKT system, whereas the 'EliminateBCsSolveSystemDirectly' did
% not. Additionally, the new Null-Space solver solves (n-m) x (n-m) system where n is the degrees of freedom and m is the
% number of constraints. Typically n is twice the number of nodes (uv solve) and m is the number of boundary conditions. The
% 'EliminateBCsSolveSystemDirectly', on the other hand, always solves an n times n system. The new solver is always faster.
% How much faster depends on the problem but one can expect it always to be at least twice as fact. Note that this only
% applies to the solution of the KKT system. If there are no boundary conditions, the same default backslash solver is used
% before.
%
% For the unsymmetrical KKT case, the solver has not changed, provided the BCs form a fat orthogonal system (This is
% typically the case if, for example, each degree of freedom is only involved in one boundary condition). If, on the other
% hand, the constraint matrix $A_{\mathrm{eq}}$ does not fulfill $A_{\mathrm{eq}} \, A_{\mathrm{eq}}^T = I_m$, a new unsymmetrical Null-Space solver is now used instead of
% the previous Augmented Lagrangian solver.
% 
% One can expect the speed-up to be noticeable in uv solves involving boundary conditions. For the uvh solve there will,
% typically, be no difference in performance, even with BCs. It is possible that if one has a large number of BCs, manually
% setting
%
%   CtrlVar.AsymmSolver="NullSpace";  
%
% might speed up the solve, as the Null Space solver solves a reduced (n-m) x (n-m) system. However, the Null-Space solver
% requires the construction of a basis for the null-space of $A_{\mathrm{eq}}$.
%
% *Release Notes* _March 2026_
% 
% An elementary additional mesh-generator is available. This one creates regular square meshes only (it is very basic
% indeed!). This can be useful when, for example, when using some simple computational domains, or to generate an initial
% mesh. All mesh refinements can then be applied to this mesh afterwards. 
%
% Example:
%
%   CtrlVar.MeshGenerator="UaSquareMesh";
%   CtrlVar.UaSquareMesh.xmin=-650e3;
%   CtrlVar.UaSquareMesh.xmax=900e3;
%   CtrlVar.UaSquareMesh.ymin=-3400e3;
%   CtrlVar.UaSquareMesh.ymax=-600e3;
%   CtrlVar.UaSquareMesh.nx=100;
%
% Note that MeshCoordinates are not used above. By specifying MeshCoordinates, all elements outside of the square are
% then automatically deactivated at the beginning of the run. An example for the use of this approach is in the Greenland
% sub-folder of Examples. 
%
% 
% 
% <<HoffSquareMeshExample.PNG>>
% 
%
%
%
% *Release Notes* _January 2026_
%
% When calculating $\dot{h}$ explicitly,  homogenized  thickness ($h$) boundary conditions are applied to the $\dot{h}$ solve. So if
% thickness is set to some prescribed value at some nodes, $\dot{h}$ at those nodes will be forced to be equal to zero. And if
% thickness links/ties are defined between nodes, those same ties will be applied to $\dot{h}$. This has no effect on $\dot{h}$
% values obtained in transient runs. But this will affect calculated/modeled $\dot{h}$ values used in an inversion when compared
% against measured $\dot{h}$ values, in combination with thickness boundary conditions.
%
%
%
% *Release Notes* _December 2025_
%
% * When using $\dot{h}$ measurements in with boundary conditions applied to h, the adjoint now uses the $\dot{h}$ calculated
% explicitly using those BCs. Previously, when evaluating the sensitivities of the $\dot{h}$-cost function term at integration
% points as part of the adjoint approach, this was not the case. 
%
% 
%
% *Release Notes* _November 2025_
%
% * The file structure of the Ua folder on github has been changed so that all key m-files are now in the sub-folder UaMain.
%
% * The mesh2d package updated to latest version.
%
% * SuiteSparse folder deleted, as it is now part of core MATLAB functionality (since R2024a)
%
% * MATLAB seems to have been busy working on their Optimisation functions, and the performance of using
%
%    CtrlVar.Inverse.MinimisationMethod="MatlabOptimisation-GradientBased";  
%
% now appears improved. This option actually stop working in Matlab 2021b and 2022a, and this may have been due to
% a bug in the optimisation toolbox. From at least 2024a onward this now works again, and based on some numerical tests,
% appears much improved. This is currently not the default option, but users might consider setting
%
%    CtrlVar.Inverse.MinimisationMethod="MatlabOptimisation-GradientBased";  
% 
% in their DefineInitialInputs.m files, to benefit from these improvements. (The old default setting using a Hessian guestimate still works as before.) 
%
% * To select the algorithm for the forward-time integration, use the string variable
%
%   CtrlVar.ForwardTimeIntegration
%
% and set to either 
%
% # "-uv-"   for a velocity solve only, i.e. no evolution of ice thickness (h). (Sometimes this is referred to as
% time-independent run, or as a diagnostic run.)
% # "-uvh-"  for an implicit velocity and thickness solve, i.e. time dependent run 
% # "-uv-h-"  for a semi-implicit velocity and thickness solve, i.e. time dependent run where the thickness is solved
% implicitly, and an outer iteration loop is used to ensure that the velocities and thickness are consistent at the end of
% the time step. This is not a recommended option, and is slower than the -uvh- option. 
%
%
% *Release Notes* _July 2025_
%
%
%
% * Enforcing post ice thickness can, as before, be done using the penalty method, by setting 
% 
%
%    CtrlVar.ThicknessPenalty=1;  
%
% This causes an additional implicit mass-balance term to be added 
%
% $$a^{\star} = a_1 (h-h_{\min}) + a_2 (h-h_{\min})^2 + a_3 (h-h_{\min})^3 $$
% 
% where 
% 
% $$h < h_{\min}$$
% 
% and where:
%
%  a_1 = CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffLin
%  a_2 = CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffQuad;
%  a_3 = CtrlVar.ThicknessPenaltyMassBalanceFeedbackCoeffCubic;
%                                                             
% Note: $a_1$ and $a_3$ need to be negative and $a_2$ positive, however, this is check internally so actually the sign on input is
% immaterial. 
%
% * By setting:
% 
%   CtrlVar.ActiveSet.ExcludeNodesOfBoundaryElements=true;
%                                                        
% nodes of boundary elements can now be excluded from the active set. By default this is false, which is the previous
% behavior.
%
% * The parameter:
%
%   CtrlVar.LevelSetMethodAutomaticallyDeactivateElementsRunStepInterval
%
% is no longer used.  Instead the interval is now identical to the general adapt-meshing interval. See further
% explanation in Ua2D_DefaultParameters.m
%
% *Release Notes* _May 2025_
%
% * The call to  
% 
%   DefineDesiredEleSize.m
% 
% now only accepts one combination of the number of input and output variables. 
% 
% * UaUtilities contains a new function: 
%
%   ElementErrorEstimator(CtrlVar,MUA,F)
%
% which provides a "Recovery-based error estimator". An example of how to use this error indicator for automated mesh
% refinement and un-refinement is provided in "1dIceStream" example within the UaExamples repository. 
%
% * Also within UaUtilities is the function 
%
%   InfluxOutfluxNodes.m
%
% which can be used to find nodes along the boundary where the velocity vectors are oriented into the domain (influx/inflow nodes) and
% out of the domain (outflux/outflow nodes). This can be used, for example, in when defining boundary condition to ensure
% that flux into the domain is some given value, e.g. zero.
%
%
% *Release Notes* _March 2025_
%
% * The uv assembly was changed slightly to make the assembly with respect to ice thickness, h, better consistent with the uvh
% assembly. This can lead to uv solve to be different if the ice density is varying greatly spatially. The basic difference
% is that the flotation condition, which involves the product rho and h, is now evaluated at integration points directly.
% Previously in the uv solve, it was evaluated by first forming the product, rho h, at nodes and then interpolating the
% product to the integration points.
%
% * Calls to gcp have been replaced with gcp('nocreate'). The implication is that a parallel pool should be defined and started
% ahead of a call to Ua. However, note that this might also depend on user settings for the parallel pool. For example if the
% user setting imply automated start of a parallel pool whenever parfor or smpd is encountered. If no parallel pool is found,
% all parallel options are turned off, and the run proceeds in non-parallel mode.
%
% * Unless the mass balance/thickness feedback is activated, MassBalance evaluation now lags behind by one time step. This was
% done to reduce then umber of calls to DefineMassBalance, and to make sure that the mass balance of the Úa fields, F0, was same as the
% F from previous time step.  Note that when the mass balance/thickness feedback is activated, the mass balance function is
% called within the assembly loop and then the mass balance will not lag behind.
%
% * The start and end times can now be specified using (new option):
%
%   CtrlVar.StartTime
%   CtrlVar.EndTime
%
% Previously this was done using (old option):
%
%   CtrlVar.time     
%   CtlrVar.TotalTime
%
% The previous option still works, but the new option is recommended. 
%
% * When not inverting for both A and C, the variable InvValues was not updated in last call, leading to a possible mismatch
% between F.A and InvValues.A. Thanks to Camilla Schelpe for identifying this.
%
% * The implementation of the ice-thickness barrier function has been simplified, and is now similar to the barrier function
% used in the level-set solver.
%
% * The MassContinuity solver, (only used when solving for h alone, and not in uv or uvh solves), now uses the active-set
% method to enforce positive thickness.
%
% The model was tested with MATLAB R2024b 
%
% *Release Notes* _September 2024_
%
% * For comparison purposes the semi-implicit solver has been updated. 
%
% This option is NOT recommended and the (fully) implicit approach continues to be the default. The semi-implicit approach is
% both slower and less accurate that the implicit one.  This option has now been updated and made available for use as
% apparently most (all?) other ice-flow models appear to be using this semi-implicit approach.
%
% This option is activated by setting:
%
%    CtrlVar.ForwardTimeIntegration="-uv-h-" ;
%
%
% * The default value of 
%
%   CtrlVar.LevelSetMethodThicknessConstraints;
%
% is now set to true (previously set to false). Therefore if one uses the active-set method and the level-set method, min
% thickness constraints will be automatically applied to all nodes downstream of calving fronts, unless of course this
% parameter is set to false by the user. 
% 
% *Release Notes* _July 2024_
%
% * The utility
%   
%   [tbx,tby,tb,eta] = CalcBasalTraction(CtrlVar,UserVar,MUA,F,options)
%
% has been updated to allow for calculation at integration points (before only calculated nodal values).
%
% * An inconsistent input parameter test when using different sliding laws for basal drag and ocean drag calculations
% corrected. Thanks to Sainan Sun for pointing this out. 
%
% * A situation where a composite variable was used inside of spmd, resulting in a warning but no errors, has been corrected.
%
% *Release Notes* _March 2024_
%
% * Positive thickness constraints --- added dynamically when using the active-set option --- are now not added to nodes already
% contained in a user-defined constraint. Hence, any user-defined thickness constraints are respected, even if this means
% that some thickness go below the minimum specified thickness. Previously, thickness constraints were added to all nodes
% with thickness less than CtrlVar.ThicMin.
% 
%
%
% *Release Notes* _February 2024_
%
% * Parallel spmd assembly option now improved and shows good scalability, although speedup always somewhat problem dependent.  
% On local workstations with 12 workers, speedup ranging from 6 to 10 seems easily obtainable, and on machines with 24 workers a
% speedup of 20 has been observed. 
%
% To switch on the smpd parallel assembly do:
%
%   CtrlVar.Parallel.uvhAssembly.spmd.isOn=true;    % uvh assembly in parallel using spmd over sub-domain (domain decomposition)  
%   CtrlVar.Parallel.uvAssembly.spmd.isOn=true;     % uv assembly in parallel using spmd over sub-domain (domain decomposition)  
%
%
% * The linear system can now also be solved using distributed arrays, although not all cases yet implemented. Turn this on/off as:
%
%   CtrlVar.Parallel.Distribute=true/false;                       % linear system is solved using distributed arrays. 
%
% Again, the speedup is problem dependent and never particularly large for low density sparse matrices as typically generated
% by Úa and other FE programs. 
%
%   CtrlVar.Parallel.isTest=false/true;             % Runs both with and without parallel approach, and prints out some information on relative performance. 
%                                                   % Good for testing if switching on the parallel options speeds things up, and by how much.
% 
% *Release Notes* _October 2023_
%
% * A (rare) case where the Cauchy direction is not a direction of descent was incorrectly updated. This has now been addressed. 
%
%
% *Release Notes* _September 2023_
%
% * Thanks to Sebastian Rosier for providing an example where the backtracking algorithm failed. This has now been corrected.
%
%
% *Release Notes* _July 2023_
%
% * uv and uvh solver now uses dog-leg search if Newton back-tracking results in small steps.
%
% * Thanks to Sainan Sun for spotting that in an adaptive mesh step, call to calving was ahead of call to geometry and
% densities. This has now been corrected.
%
% *Release Notes* _July 2022_
%
% * Call to DefineOutputs is now only done at the beginning and end of runs if the variables
%
%       CtrlVar.CreateOutputsBeginningOfRun=true;   % If true, then call DefineOutputs at the beginning of a run, that is ahead of the runstep/transient loop.
%       CtrlVar.CreateOutputsEndOfRun=true;         % If true, then call DefineOutputs at the end of a run, that is after the runstep/transient loop
%
% are set to true.
%
% * Calving options have been greatly improved and a flexible framework for implementing calving laws implemented. Although
% still considered not fully tested, this option appears to work quite well. Further details can be found in
% Ua2D_DefaultParameters.m and in the DefineCalving.m files. 
%
% * In the active-set method a minimum number of new active nodes can be specified for the active set to be updated ahead of a
% new uvh solve, e.g.:
%
%   CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints=5; 
%
%
% * UaPlots is a new plotting utility function.
%
% * Some minor bugs have been corrected, but nothing that should affect calculated solutions. 
%
% * This version requires MATLAB2021a or later. 
%
% *Release Notes* _Sept 2021_
%
%
% * Calving
%
% Calving can now be prescribed by defining a ice/ocean mask at each time step.
% The implementation is more general and faster than other similar manual
% calving options such as deactivating elements or prescribing an additional
% mass balance term.
%
% The ice-free areas are automatically melted away using a mass-balance feedback
% option implemented at the integration points. Second-order NR convergence is
% obtained even if the resulting mass balance distribution varies significantly
% spatially within an element
%
% To use this calving option set:
%
%      CtrlVar.LevelSetMethod=1;
%
% and then define the ice/ocean mask at each time step using the m-file
%
%      DefineCalving.m
%
%
%
% NOTE: Currently the only supported calving option is based on directly
% prescribing the LSF at each time step
%
%
% If you wanted, for example to get rid of all floating ice for x>500e3, do:
%
%
%
%         [UserVar,LSF,CalvingRate]=DefineCalving(UserVar,CtrlVar,MUA,F,BCs)
%           F.GF=IceSheetIceShelves(CtrlVar,MUA,F.GF);
%           OceanNodes=MUA.coordinates(:,1)>500e3 & F.GF.NodesDownstreamOfGroundingLines;
%           LSF=zeros(MUA.Nnodes,1)+ 1 ;
%           LSF(OceanNodes)=-1;
%         end
%
% Think of the LSF as a ice/ocean nodal mask, positive for ice and negative for ocean. 
%
% The the ice thickness in ice-free areas is
%
%   CtrlVar.LevelSetMinIceThickness
%
% While this could in principle be set to zero, doing so might create isolated
% islands of ice, resulting in an undetermined system. So put this to a some
% small positive value (small compared to typical ice thicknesses of interest).
% By default:
%
%    CtrlVar.LevelSetMinIceThickness=CtrlVar.ThickMin+1;
%
% where
%
%     CtrlVar.ThickMin=1;   
%
%
% Automated mesh refinement based on distance to calving from is possible in similar way as around grounding lines. For
% example to locally refine all elements with 10e3 distance units from the calving front to 1e3, use:
%
%   CtrlVar.MeshAdapt.CFrange=[10e3 1e3] ; 
%
% Remember also to set 
%
%   CtrlVar.AdaptMesh=1;
% 
% as always needed to activate the remeshing options.
%
%
% * New argument list option with F as an input
% 
% New:
%
%      DefineMassBalance(UserVar,CtrlVar,MUA,time,F)
%
% Old:
% 
%      DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
% 
% The old approach still works and is generally recommended. 
%
% 
% * Geometry and densities can now be defined in the same m-file, using 
%
%      DefineGeometryAndDensities.m 
%
% The previous approach of using separate m-files to define geometry
% (DefineGeometry.m) and densities (DefineDensities.m) still works. But if you
% have a "DefineGeometryAndDensities.m" in the run folder, only
% "DefineGeometryAndDensities.m"  is used and "DefineGeometry.m" and
% "DefineDensities.m" ignored.
%
% * x and y nodal values now a field of F, ie F.x and F.y
%
% The (x,y) coordinates of the nodes are now also accessible through F as
%
%   (F.x,F.y)
%
% These are always identical to MUA.coordinates, that is
%
%   F.x=MUA.coordinates(:,1)
%   F.y=MUA.coordinates(:,2)
%
% This was simply done for convenience, and now one can, for example, create a plot showing the upper surface, s, as a function of x as
%
%   plot(F.x.F,s,'.')
%
% * Default inverse algorithm has changed. Now the default options is a
% Hessian-based inversion.  The older-approach is still available by selecting a
% Gradient-based inversion. This is specified using the CtrlVar field:
%
%       CtrlVar.Inverse.MinimisationMethod
%
% When using this Hessian-based inversion, there is no need to specify the value of 
%
%       CtrlVar.Inverse.AdjointGradientPreMultiplier 
%
% and this field is not used.
%   
% * Ocean and wind-induced drag over floating ice can be included. This is defined in 
%
%   DefineSeaIceParameters.m
%
% and you need to switch this option on by setting
%
%   CtrlVar.IncludeMelangeModelPhysics=true;
%
% in DefineInitialInputs.m
%
% * The quadrature degree can now be directly specified using
%
%   CtrlVar.QuadratureRuleDegree=N
%
% where N is a number that can be as high as 25 (Previously the highest possible
% degree was 13). Generally, the default option should be fine, but there might
% be instances were increasing the degree could help with obtaining second-order
% convergence of the Newton-Raphson system.
%
% * The shallow-ice sheet (SIA/SSHEET) option now includes basal sliding.
%
% The SSHEET option is implemented for the Weertman sliding law only. The transient SSHEET solution is done implicitly with respect to
% the thickness using the NR method.  When using SSHEET you will, in general, need to specify boundary conditions for both the
% deformational and the basal sliding velocities. Note the SSHEET option is based on the shallow-ice approximation.
%
% *Release Notes*
% _June 2020_
%
% * Naming of some user input m-files has changed: 
%
%
%       Ua2D_InitialUserInput.m  
%
% is now named:
%
%       DefineInitialInputs.m
% 
% And
%   
%       UaOutputs.m
%
% is now:
%
%       DefineOutputs.m
%
% These changes were done so that names of all user-input files start with 'Define'
%
% Note: To systematically change the names of all your old Ua2D_InitialUserInput.m and
% UaOutputs.m files you can, for example, consider using the m-file utility file:
%
%       RenameFileRecursively.m
%
% * Naming of a few CtrlVar fields has changed. For example
%
%       CtrlVar.ATStimeStepTarget
%
% is for example no longer used. Use:
%
%       CtrlVar.ATSdtMax
%
% instead.  
%
% You can now specify a minimum selected time step in the
% automated-time-stepping (ATM) algorithm settting the field
%    
%       CtrlVar.ATSdtMin
% 
% Also
%
%   CtrlVar.UaOutputsDt 
%
% has also been replaced by
%
%   CtrlVar.DefineOutputsDt 
%
% If you define these old fields in your (new) DefineInitialInputs.m, �a will spot this and complain bitterly. 
%
% Those using Unix might want to systematically change names of some of the CtrlVar fields
% in there old input files. You might be able to use something like:
%
%   find . -name "DefineInitialInputs.m" -exec sed -i 's/CtrlVar.ATStimeStepTarget/CtrlVar.ATSdtMax/g' {} +
%   find . -name "Define*.m" -exec sed -i 's/CtrlVar.UaOutputsDt/CtrlVar.DefineOutputsDt/g' {} +
%
% * Exit criteria for the non-linear uv and uvh loops are now more flexible and allow for
% more options. See comments in 
%
%       Ua2D_DefaultParameters.m
%
% for more details.
%
% The
%
%   RunInfo
%
% variable now contains more information about the run, such as time step and number of
% non-linear iterations per runstep. 
%
% Most of the changes since last Feb 2020 have been 'under the hood'. For example solving
% a KKT type system can now be done with a much more flexible pre-elimination method than
% before.
%
% All line searches are now done using one single backtracking algorithm (before at least
% 4 different backtracking routines were used.) 
%
% The KKT system is now always solved using the primal-dual method. (Previously the
% initial iterative was made feasible and the primal method then used.) 
%
%
% *Release Notes*
% _February 2020_
%
%
% * Úa can now be called with CtrlVar as second argument, e.g 
%
%   Ua([],CtrlVar) 
%
% in which case the fields of the CtrlVar on input will be those used in the run
% even if these same fields are also defined in Ua2D_InitialUserInput.m
% 
% * The combination of local mesh refinement using the newest-vertex bisection method with
% mesh deactivation is now very flexible and allows for both reduction and extension of the
% computational domain during the run.
% 
% * Several new sliding laws are now implemented. These include Weertman, Coulomb, Budd, Tsai,
% Cornford, and  Umbi. Refer to the Úa compendium for definition of these
% different sliding laws. Inversion is also possible for all these sliding laws except for
% the Coulomb sliding.
% 
% * The semi-implicit time-stepping method now uses SUPG by default and when used the
% automated time-stepping algorithm ensures that time step fulfils the CFL condition.
% 
% * When using automated mesh refinement/coarsening in a time-dependent run, the upper
% surfaces are mapped between meshes and the thickness calculated afterwards. (Previously,
% the thickness was mapped and the surface calculated afterwards.) This may cause a
% potential violation of mass conservation. However, this also ensures smooth upper
% surfaces even when refining mesh over areas with very uneven bedrock topography.
% 
% * Inversion can be done using all implemented sliding laws except Coulomb. 
% 
%
% * u, v and h residuals now calculated in the L2 norm instead of the l2 norm as in the
% past. 
%
% *New user input file options:*
% 
%   DefineRunStopCriterion.m
% 
%   DefineFinalReturnedValueOfUserVar.m
% 
%   DefineOutsideValues
% 
% *Several new computational utilities* 
% 
%
%   LakeOrOcean3.m 
%
% is both fast and robust (Thanks to Sebastian Rosier).
% 
%   EleFlooding.m
%
% finds all nodes/elements connected to a given node fulfilling some
% additional criterion.
% 
% 
%%
