

%%
% 
% *Change Log* _October 2026_
%
% * UaFields.loadobj: objects saved with a different set of properties (2 November 2026):
%
% If a UaFields object was saved with a property that has since been removed (eg kHElement, in restart files written with the version of
% 1 November 2026), MATLAB passes it to UaFields.loadobj as a structure, which was returned unchanged, giving the warning
% MATLAB:class:mustReturnObject on loading. loadobj now creates a UaFields object and copies the fields that are still properties.
%
% * Grounding-line transition: nodal transition length, CtrlVar.kH replaced (2 November 2026):
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
% The grounding-line mask is now G=(1+tanh((h-hf)/l))/2 with a transition length l (in units of thickness, m) that is a nodal field in
% transient runs, F.lGL, evaluated once per time step from the state at its start (GroundingLineTransitionLength.m):
%
%   l = min(lmax, max(lmin, c*dt*|dh/dt|)),
%
% with dh/dt the rate of thickness change of the previous time step and dt the current time step. The transition is thus made wide
% compared with the change of thickness over the time step, which removes the limitation of the time step by the mask: the uvh Newton
% iteration converges at grounding-line Courant numbers of order one (Thule: Co_GL~1 at dt=0.5), where with lmin=1 m it stagnated from
% Co_GL~0.01. The smoothing is proportional to dt (like the error of the implicit time stepping) and vanishes as dt->0; lmax limits it
% where |dh/dt| is large far from flotation (eg near calving fronts). New options:
%
%   CtrlVar.GroundingLineTransition.lmin=1;    % (m) replaces 1/CtrlVar.kH; defining CtrlVar.kH now gives an error
%   CtrlVar.GroundingLineTransition.lmax=25;   % (m)
%   CtrlVar.GroundingLineTransition.c=0.5;     % (dimensionless); c=0 gives l=lmin everywhere (the previous behaviour with kH=1/lmin)
%
% l is kept fixed during the Newton iterations (no additional Jacobian terms), recalculated if the time step is reduced (uvh2NotConvergent.m),
% and reset on mesh changes. It is used, through IntegrationPointkH.m (kH=1/l at the integration points), by the uvh and uv assemblies
% (vectorised and compiled: the kernels take the nodal l and interpolate it), so that -uvh- and -uv-h- solve the same discrete problem, and
% by the gradient and Hessian calculations of the inverse problem. Diagnostic solves and inversions use l=lmin (no time step). The nodal
% mask of the lower surface (Calc_bs_From_hBS.m) and the water-depth mask use lmin. The element-wise kH of 1 November 2026
% (CtrlVar.kHElementAlpha, F.kHElement, ElementWisekH.m) has been removed; defining CtrlVar.kHElementAlpha gives an error.
%
% New diagnostic: GroundingLineCourantNumber.m, the grounding-line Courant number of each time step (median, 90th percentile, maximum over
% the grounding-line elements of the change of h-hf over the time step relative to its range over the element), printed at InfoLevel>=1
% and stored in RunInfo.Forward.GLCourantMedian, GLCourantP90 and GLCourantMax.
%
% Tests (Thule, fixed calving front, meshes refined to 5.8 km at the grounding line; t=0.5 to 2.5): with fixed time steps of 0.25 and 0.5 yr
% no time-step reductions (with the element-wise kH, alpha=0.02: 1 and 2; with alpha=0.05: 0 and 1); the remaining error against dt=0.01 is
% the time-discretisation error (thickness about as with alpha=0.05 at the same time step, velocities considerably better). With adaptive
% time stepping (ATSTargetIterations=25, t=0 to 5) 88 time steps and 606 Newton iterations, against 662 with alpha=0.02. Greenland (real A
% and C, dt=0.1): differences from the element-wise kH of 2 mm rms in thickness, same run time; -uvh- and -uv-h- agree to 0.6 mm rms.
% Regression tests: TestUaCompiled.m (new tests testUvhNodalTransitionLength, testUvNodalTransitionLength), TestInverseNodalTransitionLength.m
% (replaces TestInverseElementWisekH.m; adjoint gradient and Hessian consistent with a spatially varying, frozen F.lGL): 92 of 92 passed.
%
% * Time-step reduction after a failed uvh solve starts from the active set at the start of the time step (2 November 2026):
%
% When a uvh solve did not converge, uvh2NotConvergent.m reduced the time step and solved again from the velocities, thickness and
% multipliers at the start of the time step, but with the active set of thickness constraints left by the failed attempt. In a failed
% attempt with a large time step, many nodes can become (deeply) negative and are added to the active set (in addition to any predictive
% activation based on the explicit estimate for the large time step). Keeping these constraints for the smaller time step held nodes at the
% minimum thickness that would otherwise have kept their ice, and the cumulative release criterion (CtrlVar.ActiveSet.ReleaseAlpha)
% releases them only gradually. Thule (refined mesh, t=0.5 to 1.5): a run whose first step was halved from dt=1 to dt=0.5 differed from a
% run with dt=0.5 from the start by 107 m rms in thickness, with up to 2276 m of ice lost at individual nodes. uvhRootFinding.m now stores
% the boundary conditions at the start of the time step and passes them to uvh2NotConvergent.m (new optional last input), which resets the
% active set (and the related RunInfo.Forward.ActiveSetLockedNodes, ActiveSetNoFurtherReleases and ActiveSetah) to them before each new
% attempt. The two runs above now agree to round-off. This affected all runs in which a time step was reduced after a failed uvh solve.
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
% * Inverse problem with the element-wise kH (1 November 2026) [superseded by the nodal transition length of 2 November 2026, see above]:
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
% * Element-wise kH is the default, CtrlVar.kHElementAlpha=0.02 (1 November 2026) [superseded by the nodal transition length of 2 November 2026, see above]:
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
% * Compiled uv and uvh kernels with the element-wise kH (1 November 2026) [superseded by the nodal transition length of 2 November 2026, see above]:
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
% * Element-wise kH in the uv solver, computed once per time step (1 November 2026) [superseded by the nodal transition length of 2 November 2026, see above]:
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
% * uvh solver: grounding-line fallback with element-wise kH (1 November 2026) [removed again later the same day, see above]:
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
% New: element-wise kH for the mask at the integration points of the uvh assembly (uvhAssemblyIntPointImplicitSUPG_v2.m),
%
%   kH_e = min(CtrlVar.kH, 2/(CtrlVar.kHElementAlpha*range of h-hf over the element)),
%
% with the range taken from the thickness at the start of the time step, so that kH_e is fixed during the Newton solve. The derivative
% (deltaint) is evaluated with the same kH_e, so the Jacobian stays consistent (Taylor test: second order). Where the mesh resolves the
% transition, kH is unchanged; under mesh refinement the mask converges to that with the scalar kH. The nodal mask (lower surface,
% Calc_bs_From_hBS.m) keeps the scalar kH. With kHElementAlpha=0 (default) the assembly is unchanged.
%
% Grounding-line fallback (uvhRootFinding.m, on by default): if a uvh solve does not converge (and is neither an early exit on negative
% thickness nor a stall below CtrlVar.uvhStalledForceTolerance), the step is solved again, once, from the stalled iterate with
% kHElementAlpha=CtrlVar.uvhGLFallback.Alpha, before the time step is reduced. If this succeeds, the element-wise kH is kept for the
% remaining solves of this time step only; otherwise the time step is reduced as before, with the scalar kH. The compiled uvh kernel does
% not support the element-wise kH; the fallback solves use the vectorised assembly (CompiledAssemblyCheck.m). New options:
%
%   CtrlVar.kHElementAlpha=0;             % >0: always use the element-wise kH (uvh assembly)
%   CtrlVar.uvhGLFallback.isOn=true;
%   CtrlVar.uvhGLFallback.Alpha=0.05;
%
% and RunInfo.Forward.uvhGLFallbackUsed (last time step) and RunInfo.Forward.uvhGLFallbackCount (number of time steps).
%
% Thule example (fixed calving front and mesh, t=0.5 to 1.5, kH=1, -uvh-): with the fallback, all time steps up to dt=0.25 converge without
% time-step reductions (previously: reductions at all dt>=0.025); the fallback was used in 2 of 4 (dt=0.25), 5 of 10 (0.1), 6 of 20 (0.05)
% and 1 of 40 (0.025) steps. Differences from a reference (kH=1, dt=0.005) at t=1.5 are as small as those of -uv-h- with kH=1 at the same
% time step (rms thickness, dt=0.025: 1.08 m against 0.97 m; dt=0.1: 4.0 m against 4.8 m; dt=0.25: 9.2 m against 12.5 m). With the
% element-wise kH used throughout (kHElementAlpha=0.05), there is a systematic deviation of about 3.5 m rms. Greenland (-uvh-, dt=0.1),
% where all solves converge, is bit-identical (fallback not used).
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
% * uv-h outer iteration: divergence detection, and multipliers that match the active set (31 October 2026):
%
% uvhSemiImplicit.m: (1) The multipliers of the thickness constraints returned by the h solve belong to the active set of its last
% solve. If the active set was changed after that solve (the active-set loop of the h solve left at its maximum number of iterations,
% ActiveSetConverged=0, or deferred releases), they did not match BCs1, and the convergence test failed with an error in
% CalcCostFunctionNRuvh.m ("Incorrect dimensions for matrix multiplication"). They are now reset (solved for again in the next h
% solve), and such an outer iteration is not counted as converged. (2) If the residual of the outer iteration is not finite or grows by
% more than CtrlVar.uv2h.DivergenceFactor (new, default 1e3) above its smallest value (from the third outer iteration on), the outer
% iteration is stopped as diverged, and the iterate with the smallest residual is returned as not converged (RunInfo.Forward.uvhDiverged
% =true). Greenland nx150 with the inverted A and C, dt=1: the uv-h iteration diverges (rForce from 13 to 1e22, h down to -3675 m) and
% previously ended with the error above after 420 s; now stopped after 5 outer iterations. Converging runs are unchanged (bit-identical).
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
% Files: UaFields.m (new property hSuppressedThickening), UpdateSuppressedThickening.m (new), ActiveSetUpdate.m (release rule; the nodal
% rates ah of the last active-set update of each time step are kept in RunInfo.Forward.ActiveSetah, as BCs and l need not match at the
% end of a time step with deferred releases), uvhRootFinding.m and MassContinuityEquationNewtonRaphsonThicknessContraints.m (reset of
% ActiveSetah per time step), Ua2D.m (update once per accepted time step, -uvh- and -uv-h-), MapFbetweenMeshes.m (mapping, with a message
% for InfoLevelAdaptiveMeshing>=1), Ua2D_DefaultParameters.m. Tested: Greenland dt=0.01 and 0.001 (-uvh-), dt=0.01 (-uv-h-), a uniform
% mass balance of 1 m/yr at dt=0.001 (ice-free ground now thickens as at dt=0.0025), restart at t=0.25 (equal to round-off), Thule with
% deferred releases, mesh adaptation and adaptive time stepping.
%
% * ABfgPreEliminate.m: fast path for constraints that each fix one unknown (31 October 2026):
%
% The constrained system is solved as Atilde x = btilde, with Atilde=Q*A+B'B, Q=I-B'B (after scaling A). In the usual case, each row
% of B has a single entry (+-1), ie each constraint fixes one unknown (fixed values, active thickness constraints); then Q is
% diagonal, and Q*A only sets the rows of the constrained unknowns to zero. With PARDISO, and without a decomposition being reused,
% Atilde and btilde are now formed directly (elementwise multiplication with a 0/1 vector, plus a diagonal), without scaling A, and
% with the scale of the identity rows taken from the diagonal entries of the constrained rows only. Otherwise (eg ties between
% unknowns, MATLAB's decomposition, distributed arrays) the previous general path is used unchanged.
% Greenland nx150 uvh system (168,210 unknowns, 3.5 million nonzeros, 409 constraints): forming Atilde and btilde 110 -> 54 ms per
% call; x and y equal to the previous version to rounding errors (6e-16, 2e-16); general path and MATLAB solver bit-identical.
% Greenland run with -uvh- (nx150, dt=0.01, 50 steps): time steps 118.8 -> 111.0 s, results equal to 6e-12 m. -uv-h-: no measurable
% change (the h systems are smaller). Regression test: 88 of 88 pass.
%
% * CtrlVar.Parallel.BuildWorkers now true by default (31 October 2026):
%
% BuildMuaWorkers.m builds the worker data for the spmd assembly (MUA.workers) only if CtrlVar.Parallel.BuildWorkers=true. The default
% was false, and Ua2D.m set it to true ahead of each uv and uvh solve (five places), so that the worker data were not built for every
% mesh. Since 26 October 2026, UpdateMUA.m only builds them when needed (spmd on, a pool open, and the mesh, element type or pool
% changed; see MuaWorkersKey.m), so this is no longer necessary. Changes:
%   - Ua2D_DefaultParameters.m: CtrlVar.Parallel.BuildWorkers=true.
%   - Ua2D.m: the five lines setting it to true removed. (AdaptMesh.m and AdaptMesh2025.m still set it to false locally, so that no
%     worker data are built for the intermediate meshes.)
%   - BuildMuaWorkers.m: uses gcp('nocreate') and builds nothing if no pool is open (gcp would have created one, eg when called from
%     CreateMUA.m).
%   - uvhMatrixAssemblySSTREAM_SPMD.m: if there are no worker data, the assembly is done serially, as in uvMatrixAssemblySSTREAM_SPMD.m,
%     instead of failing (eg when spmd is used outside Ua2D.m without the worker data having been built).
% Thule (0.01 yr, initial mesh adaptation, level set) with spmd for uv and uvh: results bit-identical to the previous code, same
% iterations as serial; the worker data are now also built once for the initial mesh ahead of the mesh adaptation (one extra build
% per run, about 0.9 s here).
%
% * spmd versions of the uv and uvh assembly: less overhead per call (31 October 2026):
%
% Measured on Greenland nx150 (111,000 elements, 4 workers), one uv assembly with spmd took 320-334 ms (serial: 589 ms), of which
% about half was the workers' assembly and half overhead: computing the fingerprint of the fields (16-18 ms; uvh 29-32 ms), fetching
% the results (18 ms), waiting for the slowest worker (up to 27 ms), summing the sparse matrices (about 20 ms), and starting,
% synchronising and ending the spmd block (45-54 ms). Changes:
%   - SPMDFieldsFingerprint.m: the weights are kept between calls and the weighted sums are dot products (the careful treatment of
%     non-finite values only where needed): 16-18 -> 4-6 ms (uv), 29-32 -> 7-8 ms (uvh, F0 and F1).
%   - uvMatrixAssemblySSTREAM_SPMD.m, uvhMatrixAssemblySSTREAM_SPMD.m: worker 1 packs all summed results into one structure, which
%     the client fetches in a single access of the Composite (each access costs about 5 ms): 9-18 -> 6 ms. The three residual vectors
%     of the uv assembly are summed in a single spmdPlus.
% Saving about 25 ms per call (uv full 334 -> 300 ms, residual only 195 -> 176 ms). Results unchanged (spmd against serial: to
% rounding errors, also with ZeroFields and the lumped nodal penalty; Thule run with spmd: same iterations as serial). Tried and not
% kept: summing the matrix values as one dense vector on the global pattern (one spmdPlus for everything); slower, as the dense vector
% carries the whole pattern and building the sparse matrix on the client (sparse) costs 30-50 ms. With the compiled assembly, spmd is
% not needed (and not used); these changes matter only where the compiled assembly is not available.
%
% * Regression test of the compiled assembly: h tests added to the regression test (31 October 2026):
%
% UaCompiled/Tests/TestUaCompiled.m now also tests the compiled h assembly (tests testH*; 88 tests in all, about 27 s): f0 (a full
% vector) and K compared with the vectorised h assembly for level-set feedback, thickness penalty (at integration points and nodal),
% constrained nodes, mass-balance geometry feedback, all SUPG tau options, linear, quadratic and cubic elements and quadrature rules
% with 1 to 33 points, full and residual only, each time first checking that the compiled path is used; the fallback outside the
% compiled range and when dFdt is requested (identical results); one non-linear h solve with MassContinuityEquationNewtonRaphson.m
% (same iterations and solution); complex-valued input is an error. testBuildInformation covers all 13 MEX files. All 88 pass. A
% deliberate change of 1e-6 in the SUPG parameter of the vectorised h assembly only made 18 of the 23 h tests fail (all assembly
% comparisons except tau0).
%
% * Compiled h assembly (31 October 2026):
%
% CtrlVar.hMatrixAssembly.CompiledKernel=true (default false; a third switch, besides those for the uv and uvh assembly) makes
% MassContinuityEquationAssembly.m (the h assembly, used by the uv-h scheme and by MassContinuityEquation.m) assemble f0 and K with a
% compiled element-based kernel (UaCompiled/hElementKernel.m, entry points hKernelFull.m and hKernelRonly.m, MATLAB side
% hMatrixAssemblyCompiled.m), in the same way as the compiled uv and uvh assembly. The common preparation (mass-balance geometry
% feedback, level-set mask, constrained nodes) and the lumped nodal thickness penalty stay in MassContinuityEquationAssembly.m and are
% shared with the vectorised assembly. Not used if dFdt is requested (nargout=4). CompiledAssemblyCheck.m now also handles "h";
% BuildUaCompiled.m builds the h variants (hFull, hRonly, hFull_n3i6, hRonly_n3i6) as well. Complex-valued input is an error.
% Results equal to the vectorised assembly to rounding errors (at most 4e-13 relative; 2e-12 with SUPG tau1) for level-set feedback,
% thickness penalty (at integration points and nodal), constrained nodes, mass-balance geometry feedback, all tau options, linear and
% cubic elements; fallback with a message outside the compiled range and without an h assembly pattern. Greenland nx150: full
% assembly 180 -> 19.8 ms, residual only 131 -> 7.2 ms; h solves (Greenland nx150) 0.9 -> 0.35-0.45 s each, thickness equal to 1e-12 m.
% Greenland forward run (nx150, uv-h, 20 time steps), with the compiled uv, uvh and h assembly: same iterations in all time steps
% as with the vectorised h assembly, solution equal to 1e-15 relative; h assembly (332 calls) 45.3 -> 5.5 s, h solver 56.6 -> 16.9 s,
% whole run 149.1 -> 108.2 s (profiled; 313.6 s with the vectorised uv assembly with spmd and the vectorised h assembly).
%
% * h assembly with a cached assembly pattern, residual as a full vector (31 October 2026):
%
% MassContinuityEquationAssembly.m (the h assembly, used by the uv-h scheme and MassContinuityEquation.m) now builds K (and dFdt) with a
% cached assembly pattern, MUA.hAssemblyPattern, as the uv and uvh assemblies do, and returns the residual f0 as a full vector
% (accumarray) instead of a sparse one (no caller depends on it being sparse). The element matrices are accumulated as K(e,j,i), as the
% pattern expects. AssemblyPatternCache.m has a third output, hAssemblyPattern (one field, one block; switched by
% CtrlVar.MUA.AssemblyPattern.h, default true); CreateMUA.m and UpdateMUA.m create and update it like the other two, including for
% meshes read from files written before it existed. If MUA has no h pattern, the assembly creates it in each call, with a warning.
% Results: 108 of 108 outputs equal to rounding errors (largest relative difference 2.7e-16). Greenland nx150 (110,976 elements):
% full assembly 237 -> 173 ms, residual only 141 -> 128 ms. Greenland forward run (nx150, uv-h, 20 time steps): same iterations,
% solution equal to 1e-15 relative; h assembly 62.1 -> 45.3 s, whole run 163.9 -> 149.1 s.
% MassContinuityEquationNewtonRaphson.m now passes MUA.hAssemblyPattern to solveKApe.m as the structural pattern for PARDISO, as the
% uv and uvh solvers do, so that the analysis is done once per mesh whatever the pattern of the individual matrices. (Before, the
% learned pattern of the first matrix was used; for the Greenland run this also gave a single analysis.) Four h solves (Greenland
% nx150): one analysis, four factorisations, no pattern extensions, thickness bit-identical to before.
%
% * Mass balance at integration points removed (31 October 2026):
%
% CtrlVar.MassBalance.Evaluation="-int-" (evaluating the mass balance at integration points, by calling DefineMassBalance.m there) is
% no longer available, and CtrlVar.MassBalance.Evaluation is no longer defined in Ua2D_DefaultParameters.m. The mass balance is
% always evaluated at the nodes. The option was only implemented in the h assembly (MassContinuityEquationAssembly.m, used by the uv-h
% scheme and MassContinuityEquation.m), not in the uvh assembly, and was not used in any example or test. CtrlVarValidityCheck.m gives an
% informative error for "-int-"; "-node-" is accepted and has no effect. MassContinuityEquationAssembly.m no longer gathers s, b, the
% coordinates and GF at the nodes of each element (only needed for "-int-"). Results with nodal evaluation are unchanged: 108 of 108
% outputs (residual, K, dFdt; Thule and synthetic states, with level-set feedback, thickness penalty, constrained nodes and mass-balance
% geometry feedback) bit-identical. Not affected: CtrlVar.CalvingLaw.Evaluation="-int-" (calving law at integration points).
%
% * Compiled uv and uvh assembly: residual summed in the kernel (31 October 2026):
%
% The kernels (uvElementKernel.m, uvhElementKernel.m) now sum the element residuals into the nodal vectors Tint and Fext themselves,
% in parallel over the rows, using the inverse of the assembly pattern's P.iR (rPtr, rIdx, computed once per mesh in
% CompiledAssemblyPattern.m), instead of returning the element residuals for reordering and accumarray in uvMatrixAssemblyCompiled.m
% and uvhMatrixAssemblyCompiled.m. The contributions are added in the same order as accumarray: all outputs are bit-identical to
% before (112 of 112 compared, uv and uvh, all kernel variants). Greenland nx200: uv full 119 -> 100 ms, residual only 57 -> 38 ms;
% uvh full 228 -> 186 ms, residual only 90 -> 58 ms. Regression test: 65 of 65 pass. All kernel MEX files rebuilt.
%
% * Regression test of the compiled assembly: uv tests added to the regression test (30 October 2026):
%
% UaCompiled/Tests/TestUaCompiled.m now covers both assemblies (results=runtests("TestUaCompiled"), 65 tests, about 22 s; the uvh
% tests are now named testUvh*, the uv tests testUv*). The uv tests compare Ruv, Tint, Fext and Kuv of the compiled and the vectorised
% uv assembly for all sliding laws (with and without melange), linear, quadratic and cubic elements and quadrature rules with 1 to 33
% points, with and without ZeroFields, full and residual only (each time first checking that the compiled path is used, and that Ruv is
% a full vector and Kuv exactly symmetric); further the fallback outside the compiled range, one uv solve with uv.m (same Newton
% iterations and solution), and that complex-valued input is an error. All 65 pass. A deliberate change of 1e-6 in the effective
% viscosity of the vectorised uv assembly only made 23 of the 25 uv tests fail (all but the fallback and complex-input tests).
%
% * Compiled uv assembly (30 October 2026):
%
% CtrlVar.uvMatrixAssembly.CompiledKernel=true (default false) makes uvMatrixAssemblySSTREAM.m assemble Ruv, Kuv, Tint and Fext with a
% compiled element-based kernel (UaCompiled/uvElementKernel.m, entry points uvKernelFull.m and uvKernelRonly.m, MATLAB side
% uvMatrixAssemblyCompiled.m), in the same way as the compiled uvh assembly (CtrlVar.uvhMatrixAssembly.CompiledKernel). As the switch is
% in uvMatrixAssemblySSTREAM.m, it applies to all uv assemblies, including those called directly by the inverse code (adjoint). The
% element matrices are summed into the nonzeros and symmetrised inside the MEX file. Not used within spmd workers
% (uvMatrixAssemblySSTREAM_SPMD.m switches it off there). Complex-valued inputs are errors.
% Greenland nx200 (197,204 elements): full assembly 1,494 -> 127 ms, residual only 806 -> 56 ms. Results equal to the vectorised
% assembly to rounding errors for all sliding laws, with and without melange, linear, quadratic and cubic elements, with and without
% ZeroFields; Thule run (0.01 yr, mesh adaptation): same uv and uvh Newton iterations, solution equal to 1e-15 relative.
% Greenland inversion (nx100, A, B, C, 2 iterations) with the compiled uv assembly against the vectorised one (with spmd): identical
% cost function history (J, I, R), estimates equal to 1e-16, gradients (adjoint) to 3e-14; wall time 42.3 -> 27.6 s. Greenland forward
% run (nx150, 110,976 elements, uv-h scheme, 20 time steps): same iterations in all time steps, solution equal to 1e-15 relative; initial
% uv solve 24.3 -> 10.0 s, whole run 307 -> 168 s.
% Threads (OMP_NUM_THREADS), Greenland nx200, full / residual-only uv assembly: 1: 1,065 / 524 ms; 4: 361 / 172; 8: 204 / 108;
% 16: 135 / 74; 32 (default): 112 / 55 ms (vectorised, serial: 1,494 / 806 ms).
% Profiles: inversion nx100: uv assembly 16.3 s (spmd) -> 3.4 s, whole run 45.2 -> 29.2 s; the largest item is then
% WriteAdjointRestartFile (10 s). Forward run nx150 (uv-h, 20 time steps): uv assembly (593 calls) 166.5 s (spmd, of which 133.6 s
% spmd machinery) -> 29.8 s, whole run 313.6 -> 168.7 s; the largest item is then the (not compiled) h assembly
% MassContinuityEquationAssembly (59 s). Of the 29.8 s, about 8 s are the summation of the residual on the MATLAB side
% (uvMatrixAssemblyCompiled.m), which could be done in the kernel.
%
%   - CompiledAssemblyCheck.m (replaces uvhCompiledKernelCheck.m) decides for both assemblies (third input "uv" or "uvh") whether
%     the compiled assembly can be used; messages, up-to-date check and switching off after a failure are separate for uv and uvh.
%   - CompiledAssemblyPattern.m (replaces uvhCSCPattern.m) gives the per-mesh data for both, with one cache entry per assembly.
%   - BuildUaCompiled.m builds the uv variants (uvFull, uvRonly, uvFull_n3i6, uvRonly_n3i6) as well; UaCompiledBuildInfo.m records
%     the source files per kernel.
%   - Ua2D_DefaultParameters.m: CtrlVar.uvMatrixAssembly.CompiledKernel, CtrlVar.uvMatrixAssembly.Ronly and
%     CtrlVar.uvAssembly.ZeroFields now have defaults (the latter two were always set by the callers).
%   - uvMatrixAssemblySSTREAM.m: CtrlVar.uvAssembly.SkipSymmetrisation and the special case for complex-step differentiation removed;
%     the tangent matrix is always symmetrised.
%
% * uvMatrixAssemblySSTREAM.m: the uv residual now full (30 October 2026):
%
% Tint and Fext, and so Ruv (and fext0 in SSTREAM2dNR2.m), are now full vectors (accumarray), as in uvhMatrixAssembly.m; they were
% sparse vectors with almost all entries non-zero. No caller depends on them being sparse. The values change at the level of
% rounding errors only (order of summation). Greenland nx200 (197,204 elements): full assembly 1,639 -> 1,457 ms, residual only
% 914 -> 751 ms; Kuv bit-identical. Thule run (0.01 yr, with mesh adaptation): same uv and uvh Newton iterations in all time steps,
% solution equal to 1e-15 relative. Greenland inversion (nx100, uv assembly with spmd): identical cost function history, estimates
% of A, B and C bit-identical.
%
% * Regression test of the compiled uvh assembly: UaSource/UaCompiled/Tests/TestUaCompiled.m (26 October 2026):
%
% A matlab.unittest class, run with results=runtests("TestUaCompiled"). It builds a small synthetic state itself (grounded and
% floating ice, thickness near ThickMin, a level-set field, melange fields, non-converged velocities), and checks, for all sliding
% laws (with and without melange), SUPG tau options, thickness penalty functions, the thickness barrier, active-set constraints,
% level-set and mass-balance geometry feedback, ResetThicknessToMinThickness, Enforce_bAboveB, linear, quadratic and cubic elements
% and quadrature rules with 1 to 33 points, that the compiled path is actually used and that residual and Jacobian agree with the
% vectorised assembly (relative tolerance 1e-11, full and residual only); further OnlyCalcBasalDragAndEffectiveViscosity, the
% fallback outside the compiled range, one implicit uvh time step (same Newton iterations and solution), and that the MEX files are
% up to date. 40 tests, about 25 s, all passing. A deliberate change of 1e-6 in the SUPG parameter of the vectorised assembly only
% made 35 of them fail.
%
% * Compiled uvh assembly: outdated or failing MEX files fall back to the vectorised assembly (26 October 2026):
%
%   - BuildUaCompiled.m records, for each MEX file, the MATLAB release and a checksum of each source file compiled into it, in
%     UaCompiled/UaCompiledBuildInfo_<mexext>.mat (new UaCompiledBuildInfo.m; checksums with ArrayChecksum.m, of the file contents
%     without carriage returns, so that other line endings do not count as a change). uvhCompiledKernelCheck.m compares these once per
%     session and variant (about 0.1 s): if the MEX files were built with another release, or a source file changed since, the
%     vectorised assembly is used, with a message naming the files. Without build information, a note says that this could not be
%     checked.
%   - A failure of the MEX files themselves (eg built for another MATLAB release or platform: "Invalid MEX-file", or an outdated
%     interface) is reported as UaCompiled:KernelFailed by uvhMatrixAssemblyCompiled.m; uvhMatrixAssembly.m catches only this error,
%     reports it once, switches the compiled kernel off for the rest of the session, and uses the vectorised assembly. Other errors
%     (eg NaN in the matrix) are passed on as before.
% Tests: unchanged sources: compiled path used; other line endings: still up to date; a changed source file and another release:
% detected, vectorised assembly used; kernel MEX file replaced by another MEX file and by an invalid file: error reported once,
% vectorised assembly used, results identical to the vectorised assembly.
%
% * Inverse run (Greenland, nx100, 2 iterations of an A, B, C inversion) with the rewritten BasalDrag.m against the code before the
% rewrite (commit 0126b5d): identical cost function history (J, I, R), final estimates of A, B and C equal to rounding errors
% (1e-16 relative), gradients to 2e-14.
%
% * SUPGtau.m and ThicknessPenaltyMassBalanceFeedback.m usable by MATLAB Coder; the compiled kernel uses them directly (26 October 2026):
%
% The compiled uvh assembly (UaSource/UaCompiled) used modified copies of these functions (SUPGtauFlex.m,
% ThicknessPenaltyMassBalanceFeedbackFlex.m), which are now removed; the kernel calls the originals, as it does BasalDrag.m.
%   SUPGtau.m: the option can be given as the name (CtrlVar.uvh.SUPG.tau) or as an integer code (new SUPGtauCode.m, also used by
%     MakeKernelParamsFlex.m, so that strings are not compared at every integration point of the compiled kernel); the third output
%     K is defined for all options (NaN except for "tau1"); plotting excluded from compiled code (coder.target('MATLAB')). Bug fix:
%     the plotting (CtrlVar.doplots and CtrlVar.PlotSUPGparameter) tested ~isempty(MUA), but MUA is not an input, so it failed.
%   ThicknessPenaltyMassBalanceFeedback.m: the width warning (with its persistent variable) and the diagnostic output excluded from
%     compiled code. Bug fix: with InfoLevelThickMin>=10, the example curves were plotted for all penalty functions, but only exist
%     for "softplus", so the diagnostic output failed for "polynomial" and "exponential".
% In MATLAB the results are bit-identical to before (all options, name or code). Tests: compiled against vectorised for the five
% tau options, the three penalty functions and the seven sliding laws: agreement to rounding errors; Thule example runs.
% Note: for "tau1", K=coth(ECN)-1/ECN loses accuracy by cancellation for small ECN (the Taylor expansion is only used below 1e-6).
%
% * Worker data for the spmd assembly (MUA.workers): reliable and cheap check if they have to be built (26 October 2026):
%
% UpdateMUA.m is called twice in each run step. With spmd assembly, it called BuildMuaWorkers.m each time, which checked the existing
% worker data by transferring the worker data of the first worker to the client (MUAworkers{1}.Nnodes) and comparing the element
% partitions on the workers in an spmd block: about 45 ms per call for a mesh of 110,000 elements (Greenland nx150), and it did not
% detect a change of the mesh that leaves the numbers of nodes and elements unchanged (eg nodes moved, or elements renumbered), as the
% partitions only depend on these numbers; the workers would then have kept the old mesh.
% Now:
%   - UpdateMUA.m keeps a key of what the worker data were built for (MUA.workersKey, see the new MuaWorkersKey.m: the key of the mesh
%     (MeshKey.m), element type, quadrature rule, number of workers), and checks: the key, numel(MUA.workers) equal to the number of
%     workers of the pool, and all(exist(MUA.workers)), which is false for all workers if the pool has been deleted or replaced by a
%     new one, even of the same size (the answer to the question in the comment there). Only if one of these does not hold are the
%     worker data built. This takes well under a millisecond, and does not communicate with the workers.
%   - BuildMuaWorkers.m no longer checks, it builds (unless CtrlVar.Parallel.BuildWorkers=false); the callers decide when.
%   - CreateMUA.m sets MUA.workersKey when it builds the worker data.
%   - UpdateMUA.m: no second build/check when the mesh has changed, and the element areas are only recalculated if needed.
% Tests: with a 4-thread pool, unchanged mesh: no rebuild, UpdateMUA 2.3 ms (Thule test mesh); nodes moved, elements renumbered, new
% pool of the same size, pool of another size: rebuilt, and the spmd assembly then agrees with the serial one (to rounding errors);
% BuildWorkers=false respected. Greenland nx150, 1 yr: UpdateMUA 8.5 s -> 4.4 s (67 calls), BuildMuaWorkers 67 -> 4 calls, results
% bit-identical. Thule with mesh adaptation and spmd for uvh and uv: worker data built twice (once for each mesh used in the solves),
% same time steps and Newton iterations as the serial reference, results identical to rounding errors.
%
% * BasalDrag2.m renamed to BasalDrag.m, and unused files removed (26 October 2026):
%
% The rewritten basal drag (see the entries below) replaces the earlier version under its original name: BasalDrag2.m ->
% BasalDrag.m, BasalDrag2dC.m -> BasalDragdC.m; all calls and references (UaMain, UaCompiled) changed accordingly. Note that the
% call of BasalDrag.m has changed: BasalDrag(CtrlVar,Derivatives,He,delta,h,H,rho,rhow,g,ub,vb,C,m,q,muk,V0,uo,vo,Co,mo,ua,va,Ca,ma).
% Removed: BasalDragOld.m, BasalDragPreUa2026b.m, rpCWN0.m, rCWN0.m, rCWV0.m (no longer used), uvhAssemblyIntPointImplicitSUPG.m
% (replaced by uvhAssemblyIntPointImplicitSUPG_v2.m), uvMatrixAssemblySSTREAM_Parfeval.m and
% uvMatrixAssemblySSTREAMpartitionTriplets.m (no longer used, see the entry of 17 October). Also removed: the setting of
% CtrlVar.BasalDrag.CalculateDerivatives in uvhMatrixAssembly.m, uvhAssemblyIntPointImplicitSUPG_v2.m and uvMatrixAssemblySSTREAM.m,
% and CtrlVar.Inverse.dFuvdClambda in Ua2D_DefaultParameters.m, both no longer read. MEX files in UaCompiled rebuilt.
% Tests: outputs of all callers bit-identical to before the rename; compiled against vectorised for all sliding laws as before;
% Thule and RadialIceCap examples run.
%
% * Compiled uvh assembly uses BasalDrag2.m directly; examples use the new sliding-law aliases (26 October 2026):
%
% UaSource/UaCompiled: the kernel (uvhElementKernel.m) calls BasalDrag2.m from UaMain, which MATLAB Coder compiles directly, with
% an explicit derivative switch (0 if Ronly, otherwise 2). BasalDragFlex.m (the modified copy of BasalDrag.m) removed.
% MakeKernelParamsFlex.m: the sliding-law code from SlidingLawCode.m instead of its own table (which still had "W" and "C"), and
% the check for the inverse branch (CtrlVar.Inverse.dFuvdClambda) removed. In the kernel, q, muk and V0 start as NaN and the
% melange variables as 0 instead of [] (fixed types). All sliding laws and aliases are now supported by the compiled kernel.
% MEX files rebuilt. Tests: compiled against vectorised for all 7 sliding laws, with and without melange, full and residual only:
% agreement to rounding errors; Thule full run with the compiled kernel: same time steps and Newton iterations as the vectorised
% reference, max |dh| 5e-12 m; Greenland nx200, one assembly: 239 ms (before: 244 ms), residual only 90 ms (93 ms).
%
% UaExamples: "W" -> "-W-" (MismipPlus, RadialIceCap), and the remaining comparisons with sliding-law names in the examples
% (Calving, Greenland, Inverse, PIG-TWG) use SlidingLawCode.m; the Inverse example tested for "rCW-V0" (capital V).
%
% * Sliding-law names and aliases: one definition, SlidingLawCode.m (26 October 2026):
%
% The names of the sliding laws were compared against string lists in many places, with inconsistent aliases. All now use
% SlidingLawCode.m: CtrlVarValidityCheck.m (which calls it first, so that the former aliases "W", "B", "C" give an informative
% error), Ua2D_DefaultParameters.m (CtrlVar.MustBe.SlidingLaw updated to the accepted names, and the description of the laws),
% DefineInputsForInverseRun.m, GetSlipperyDistribution.m, TestSlipperinessInputValues.m, dFuvdC.m, dIdCq.m, FBC.m, FCuv.m, Fqq.m,
% WeertmanSecondOrderDerivatives.m, TestInverseHessianEndToEnd.m, and the template DefineInitialInputs.m ("W" -> "Weertman").
% Inconsistencies found and corrected on the way:
%   - Joughin: "rCW-V0" (capital V) in the MustBe list and in TestSlipperinessInputValues.m, but "rCW-v0" in BasalDrag.m, so
%     the alias did not work (and V0 was set to []). Now "rCW-v0" throughout.
%   - GetSlipperyDistribution.m required muk (5 outputs) for "W" (Weertman) and "W-N0" (Budd), but not for "Weertman". Now
%     q is required for Budd and muk for Coulomb, Tsai, Cornford and Umbi, whatever the alias.
%   - dFuvdC.m and dIdCq.m used q only for "Budd", not for its alias "W-N0"; FCuv.m accepted only "Weertman", not "W".
%   - WeertmanSecondOrderDerivatives.m: its error message for other sliding laws was a string array, so that error() itself
%     failed ("badMessageArgument"), eg in FBB.m.
%   - DefineInitialInputs.m described "rCW-N0" as Cornford; it is Umbi (Cornford is "rpCW-N0").
% Tests: all accepted names pass CtrlVarValidityCheck, the former aliases give the informative error; the outputs of the BasalDrag2
% callers are unchanged for all sliding laws; Thule example runs.
%
% * BasalDrag2.m: rewritten basal drag (26 October 2026):
%
% BasalDrag.m had grown by adding one sliding law after another. New, replacing it:
%   SlidingLawCode.m : the one place where the names and aliases of the sliding laws are defined (integer code 1-7). New aliases
%                      "-W-" (Weertman), "-B-" (Budd), "-C-" (Coulomb); "W", "B" and "C" are no longer accepted (error with a
%                      message naming the new alias).
%   SlidingLawBeta.m : all sliding laws are isotropic, taub = beta (ub,vb); for each law the drag coefficient beta(U,h,C) and its
%                      derivatives with respect to U, h and C, as explicit expressions. Replaces the local functions of BasalDrag.m
%                      and the generated rpCWN0.m, rCWN0.m and rCWV0.m (Cornford, Umbi, Joughin).
%   BasalDrag2.m     : drag and its derivatives, assembled from beta in the same way for all laws, plus the melange terms. Which
%                      derivatives are calculated is given explicitly (Derivatives=0,1,2: none, u and v, u v and h), not by the number
%                      of output arguments or CtrlVar.BasalDrag.CalculateDerivatives. All outputs are always defined. MUA and B are
%                      no longer inputs.
%   BasalDrag2dC.m   : the derivative with respect to C (inverse problem), instead of returning it in taubx when
%                      CtrlVar.Inverse.dFuvdClambda=true.
% Written so that MATLAB Coder can use them directly (tested: compiled results identical to MATLAB). NRosier (could not run) removed.
%
% Callers changed to BasalDrag2/BasalDrag2dC: uvhAssemblyIntPointImplicitSUPG_v2.m, uvMatrixAssemblySSTREAM.m, dFuvdBGeneral.m,
% dIdBqGeneral.m, FBB.m, dFuvdC.m, dIdCq.m, CalcBasalTraction.m (the setting and resetting of CtrlVar.Inverse.dFuvdClambda, and of
% CtrlVar.BasalDrag.CalculateDerivatives in the inverse routines, removed). BasalDrag.m itself is not yet removed.
%
% Tests: BasalDrag2 against BasalDrag for all 7 sliding laws, m=1 and 3, with and without melange: agreement to 1.6e-15 (relative);
% BasalDrag2dC against the dFuvdClambda branch: 1.2e-15; derivatives against finite differences; outputs of the 8 changed callers
% before and after, for all sliding laws: agreement to rounding errors (at most 1.3e-14). Thule example, full run with the
% vectorised assembly: same time steps and Newton iterations, solution identical to rounding errors (max |dh| 1e-11 m).
%
% * Compiled uvh assembly: less data movement into and out of the MEX files (25 October 2026, later):
%
% Timing experiments showed that the arithmetic in the compiled kernel scales well with the number of threads, and that about half of the
% kernel time was spent moving data into and out of the MEX file. Changes (UaSource/UaCompiled, and the call in uvhMatrixAssembly.m):
%   - MUA.Deriv and MUA.DetJ are passed as separate arrays instead of as fields of a structure: structure fields are copied into a MEX
%     file, arrays are not (54+9 MB at nx200 in the Greenland example).
%   - basal drag, effective viscosity and floating mask at the integration points are only calculated and returned if requested
%     (uvhMatrixAssembly.m called with more than 4 outputs, eg by CalcBasalTraction.m).
%   - outputs that are completely overwritten are not zero-filled (coder.nullcopy).
%   - the element-matrix values are summed into the nonzeros of the sparse matrix within the kernel MEX file, so they are no longer
%     returned (122 MB at nx200); GatherSumCSC.m and GatherSumCSC_mex removed.
% Results are bit-identical to before. Greenland nx200, one assembly: residual and Jacobian 386 -> 244 ms, residual only 159 -> 93 ms;
% ratio serial vectorised : spmd vectorised (4 threads) : compiled = 16.0 : 6.9 : 1 (residual and Jacobian), 15.4 : 5.1 : 1 (residual
% only). Thule example, full run: uvh assembly 24 s -> 17 s (vectorised: 91 s). Greenland nx200, full run (0.1 yr) with the compiled
% uvh assembly and spmd for the uv assembly: 133 s (spmd for both: 189 s). The MEX files must be rebuilt (BuildUaCompiled.m), as the
% interface of the kernel has changed.
%
% * uvh assembly: optional compiled element-based assembly (25 October 2026):
%
% New option CtrlVar.uvhMatrixAssembly.CompiledKernel (default false, ie no change unless set). If true, and the options, the element
% type and the cached assembly pattern allow it, uvhMatrixAssembly.m calls uvhMatrixAssemblyCompiled.m (new folder UaSource/UaCompiled,
% see README.md there) instead of the vectorised assembly. Otherwise, or if UaCompiled is not on the path, the vectorised assembly is used,
% unchanged; uvhCompiledKernelCheck.m decides this for each call and prints the reason once if the compiled assembly cannot be used.
%
% The compiled assembly consists of:
%   - a kernel (uvhElementKernel.m) that loops over elements, with the loop over integration points inside, generated with MATLAB Coder
%     and parallelised over elements with OpenMP. It gathers the nodal values of each element itself (from the nodal vectors and the
%     connectivity), and contains the terms of uvhAssemblyIntPointImplicitSUPG_v2.m and uvhNodalLoopSSTREAM_v2.m. Options (sliding law,
%     level-set feedback, penalty, barrier, melange, SUPG tau, ...) and all parameters are run-time inputs (MakeKernelParamsFlex.m), so no
%     rebuild is needed when they change. MEX files: one for any element type with at most 10 nodes and 37 integration points, and a
%     faster one for linear elements with 6 integration points; each for the residual and Jacobian, and for the residual only.
%   - the element-matrix values are summed directly into the nonzeros of the sparse matrix (GatherSumCSC, parallel over the nonzeros,
%     using an inverse of MUA.uvhAssemblyPattern.map computed once per mesh in uvhCSCPattern.m), and the sparse matrix is created directly
%     from its CSC arrays (SparseFromCSC.c, hand-written C MEX), replacing accumarray and sparse. Entries that sum to exactly zero are kept
%     as stored zeros (sparse drops them); the values are otherwise identical.
%   The kernel calls modified copies of BasalDrag.m, SUPGtau.m and ThicknessPenaltyMassBalanceFeedback.m (in UaCompiled), with integer
%   codes instead of strings for the sliding law and the SUPG tau option, and with plotting excluded from the compiled code. Not compiled:
%   the inverse-run branch of BasalDrag.m and the sliding-law aliases "B" and "-C-" (not accepted by the forward branch of BasalDrag.m).
%
%   uvhMatrixAssembly.m: the switch, after the common preparation and before the vectorised assembly (which is unchanged).
%   uvhMatrixAssemblySSTREAM_SPMD.m: the compiled assembly is not used within spmd workers (its OpenMP threads would compete with them).
%   Ua2D_DefaultParameters.m: CtrlVar.uvhMatrixAssembly.CompiledKernel=false.
%   UaSource/UaCompiled: the compiled assembly, the MEX files and the build script BuildUaCompiled.m.
%
% Tests: with CompiledKernel=false the results are bit-identical to before. With CompiledKernel=true residual and Jacobian agree with the
% vectorised assembly to rounding errors for all supported options, for linear, quadratic and cubic elements and quadrature rules with 1 to
% 33 points. Thule example, full run (0.1 yr, 51 time steps, mesh adaptation, level set, active set): identical time steps and Newton
% iterations, solution identical to rounding errors (max |dh| 9e-12 m); uvh assembly 91 s -> 24 s, total 209 s -> 152 s (profiled).
% Greenland example, nx200 (197,204 elements), one assembly, ratio serial vectorised : spmd vectorised (4 threads) : compiled =
% 10.1 : 4.2 : 1 for residual and Jacobian (3894 / 1628 / 386 ms), and 8.9 : 3.1 : 1 for the residual only. For larger meshes, the
% fastest combination was the compiled uvh assembly with spmd for the uv assembly (CtrlVar.Parallel.uvAssembly.spmd.isOn=true).
%
% * Active set: option to defer releases to the next time step (25 October 2026):
%
% New option CtrlVar.ActiveSet.DeferReleases (default false, ie no change unless set). If, after a converged uvh solve, the active-set
% update only releases nodes (no activations), the solution of that solve (with these nodes still at ThickMin) is accepted, and the
% active-set loop is left without a further solve; the nodes are released at the start of the next time step. Activations are always
% followed by a further solve, so h>=ThickMin is still enforced in each time step; only the sign condition on the multipliers lags by one
% time step. CtrlVar.ActiveSet.DeferReleasesMaxThickening (default inf): releases are only deferred if, for all released nodes, the
% estimated thickening within the time step (-ah*dt) is less than this times ThickMin. Changes:
%
%   uvhRootFinding.m: the decision, after the call to ActiveSetUpdate. The deferred nodes are kept, with the topology key of the mesh, in
%   RunInfo.Forward.ActiveSetDeferredReleases, and their number and estimated thickening are printed.
%   ActiveSetInitialisation.m: deferred nodes are excluded from the predictive activation (and from the initial active set if that is
%   rebuilt), as they are at ThickMin and would otherwise be activated again right away. Ignored if the mesh has changed (new local
%   function DeferredReleases).
%   ActiveSetUpdate.m: the estimated thickening of the released nodes is stored in RunInfo.Forward.ActiveSetReleasedThickening
%   (diagnostic only).
%   Ua2D_DefaultParameters.m: the two new options.
%
% Motivation: in the Thule example (50 km, 51 time steps) 40 of the 43 time steps that needed more than one uvh solve did so only
% because nodes were released, and in 9 time steps some of the released nodes fell below ThickMin again (and were locked), requiring a
% third solve. With DeferReleases=true: uvh solves 106 -> 66 (time steps with a single solve 8 -> 41), uvh Newton iterations 404 -> 312,
% locking events 9 -> 2, total run time 181 s -> 153 s; time steps unchanged. Effect on the solution after 0.1 yr: relative change in
% ice volume -6e-6, RMS change in h 0.8 m, max 21 m. The changes are concentrated at the calving front (149 of the 161 nodes with
% |dh|>0.1 m are floating, within about 40 km of nodes at ThickMin), where h changes by about 20 m per time step, and are everywhere
% smaller than the change in h over one time step at that node. The mass removed by the thickness constraints increased by 6%.
% The estimated thickening of the deferred nodes was mostly just above the release threshold (median about 1.1e-3*ThickMin), but up
% to 12*ThickMin for single nodes.
% With DeferReleasesMaxThickening=1 (deferral only if no released node is estimated to thicken by more than ThickMin): 17 deferrals
% (33 nodes), uvh solves 106 -> 89, run time 181 s -> 171 s, locking events unchanged (9); relative change in ice volume -2e-7, RMS
% change in h 2 cm, max 0.6 m, mass removed by the constraints +0.5%. Thus the releases with large estimated thickening account for most
% of both the saving and the change in the solution.
% Greenland example (nx100, 25,095 nodes, 0.1 yr, 9 time steps, no calving-front evolution): 3 deferrals (52 nodes), all with small
% estimated thickening (at most 0.04*ThickMin), so DeferReleasesMaxThickening=1 and inf give identical results; uvh solves 15 -> 12,
% Newton iterations 53 -> 49 (about 5%); relative change in ice volume 1e-8, RMS change in h 0.7 mm, max 8 cm (at a node with h=1.7 m).
%
% * Level-set Newton-Raphson: fewer assemblies and iterations (25 October 2026):
%
% In each iteration of LevelSetEquationNewtonRaphson.m the level-set equation was assembled four times, each time with the matrix: once
% for the Newton system, and three times in the line search (CalcCostFunctionLevelSetEquation.m at gamma=0, at gamma=1, and again at
% the gamma returned by BackTracking.m). Changes:
%
%   LevelSetEquationAssemblyNR2.m, LevelSetEquationAssemblyNR2consistent.m: only the residual is assembled if only the residual is
%   requested. LevelSetEquationAssemblyNR2consistentScalar.m already skipped the matrix for nargout==2, but the two wrappers always
%   requested it. LevelSetEquationAssemblyNR2consistentScalar.m: K=(K+K')/2 (no advection term) only if K is requested.
%   CalcCostFunctionLevelSetEquation.m: optional input R, the residual at gamma=0 if already assembled.
%   LevelSetEquationNewtonRaphson.m: the residual assembled for the Newton system is used for the cost function at gamma=0, and the
%   cost function at gamma=1 is not evaluated again if BackTracking.m returns gamma=1.
%
% Each iteration now needs one assembly with the matrix and, if gamma=1, one of the residual only. The results are unchanged (identical
% iteration output). Thule example (50 km, 51 time steps, 102 level-set iterations): level-set assemblies 408 -> 204 (102 of them residual
% only), time in LevelSetEquation.m 45 s -> 21 s, total run time 206 s -> 181 s.
%
% LevelSetEquationNewtonRaphson.m: the hard-wired rule forcing a further iteration if rRatio<0.01 and rWork>1e-20 has been removed. It
% forced a second iteration in every level-set solve of the Thule example, although the first iteration reached r/r0 of about 1e-13
% (the second took rForce to about 1e-26). Convergence is now decided by the LSFAcceptable/LSFDesired tolerances alone. Thule example:
% level-set iterations 102 -> 77 (25 solves converge after one iteration; in the other 26 rWork is still above
% LSFAcceptableWorkAndForceTolerances(1)=1e-10 after the first), time in LevelSetEquation.m 21 s -> 16 s. The time steps and the uvh
% iteration output are unchanged.
%
% * New LinesolveScalingTest.m: performance of the linear solvers (MATLAB and PARDISO) as a function of the number of cores (17 October 2026):
%
% For a restart file, the uvh system (constraints eliminated as in ABfgPreEliminate.m) or the uv system (nodal constraints removed,
% Cholesky) is built once and saved; each measurement is then a separate MATLAB process restricted to p cores by processor affinity
% (start /affinity on Windows, taskset on Linux), as PARDISO with the TBB layer ignores MKL_Set_Num_Threads, so that its number of
% threads can not be limited from within MATLAB. Measured: PARDISO's phases separately (analysis, factorisation, substitution), with a
% new analysis for each solve and with the analysis reused; MATLAB's backslash, and its factorisation (lu/chol) and substitution; both
% solvers as called in Ua (solveKApe/solveKApeSymmetric). Model T(p)=a+b/p+dp, recommendation, plot.
%
% Results (Greenland nx150 uvh system, 168,210 unknowns, 16 cores): PARDISO with the analysis reused 1.05 s on 1 core, 0.21 s on 16
% cores (factorisation 6.6 times faster); the analysis (about 1.4 s) does not parallelise; MATLAB's backslash about 2.0-2.3 s with any
% number of cores. Unrestricted: PARDISO 1.3-1.4 times faster than MATLAB with all three phases, about 10 times with the analysis
% reused, 6.5 times as called in Ua; the analysis is 84% of a solve with all three phases. Most of PARDISO's advantage in Ua therefore
% comes from reusing the analysis, which the fixed (structural) sparsity pattern of each mesh allows.
%
% * Parallel uv assembly revised (17 October 2026), along the lines of the uvh assembly (16 October 2026):
%
% uvMatrixAssemblySSTREAM_SPMD.m was correct, but its overheads made it hardly faster than the serial assembly (Greenland nx150, 4 thread
% workers: 0.56 s instead of 0.62 s for the residual and the Jacobian, 0.41 s instead of 0.32 s for the residual only, ie slower): in each
% call CtrlVar, F and BCs were sent to all workers (0.17 s), and each worker calculated the uv assembly pattern of its sub-domain (and the
% uvh pattern) anew (about 0.23 s), as the worker data contained no uv pattern. Changes:
%
%   uvMatrixAssemblySSTREAM_SPMD.m: F, CtrlVar and BCs are kept on the workers between calls. Only ub and vb, and the CtrlVar fields that are
%   switched between calls (uvAssembly, uvMatrixAssembly, NRuvIncomplete) are sent in each call; F is sent again if the fingerprint of the
%   mesh and the other fields changes (eg AGlen or C between the forward solves of an inversion, or after remeshing), CtrlVar and BCs if
%   they change (isequaln). One spmd block for the assembly and the summation (previously two). Tint and Fext are summed over the workers
%   as well (previously returned empty). The number of workers is CtrlVar.Parallel.uvAssembly.spmd.nWorkers (or the pool size), and the
%   number of sub-domains in MUA.workers is used if they differ (previously CtrlVar.Parallel.uvhAssembly.spmd.nWorkers was set instead).
%   Without a pool or without worker data, the serial assembly is used. Optional fifth output Timing, as in the uvh assembly.
%   BuildMuaWorkers.m: the worker data contain the uv assembly pattern of the sub-domain as well as the uvh pattern.
%   New SPMDFieldsFingerprint.m: the fingerprint of the inputs kept on the workers, previously local functions of
%   uvhMatrixAssemblySSTREAM_SPMD.m, now shared by the uv and the uvh assembly (the key of the uvh assembly is unchanged).
%   New uvAssemblyScalingTest.m: the uv counterpart of uvhAssemblyScalingTest.m (restart files).
%
% Results (nx150, 56,070 nodes, 16 cores, thread pool): residual and Jacobian 0.34 s instead of 0.61 s (4 or 8 workers), residual only
% 0.21 s instead of 0.32 s; of the 0.34 s, 0.24 s are the assembly on the workers, 0.016 s the fingerprint and 0.003 s sending ub and vb.
% Results identical to the serial assembly (K to 2e-16, R, Tint and Fext to 3e-15 relative to ||Fext||). Changes of ub and vb, of
% another field (C), of CtrlVar (SpeedZero) and of the mesh (reordered elements) are handled correctly. uv solve (perturbed nx150 state):
% same 14 iterations, 19.3 s instead of 23.9 s, velocities identical to 2e-12 m/a. uvAssemblyScalingTest recommends 4 thread workers
% (1.7 times faster assembly), as uvhAssemblyScalingTest for the uvh assembly. The tests of the parallel uvh assembly give the same
% results as before.
%
% * New TestuvJacobian.m: Taylor test of the Jacobian of the uv system (17 October 2026):
%
% The uv counterpart of TestuvhJacobian.m, with the same method, options and output: for a state x=[ub;vb], the remainder
% ||R(x+eps d)-R(x)-eps K d|| is calculated for decades of eps, with R and K from uvMatrixAssembly.m (as in SSTREAM2dNR2.m), for the
% Newton direction (calculated at the state, with the boundary conditions) and smooth random directions of u, v and both. Called with
% the name of a restart file, or with (CtrlVar,MUA,F,BCs); other states are tested by changing F.ub and F.vb. In addition, the
% symmetry of K is reported. Option ModifyJacobian (a function applied to K before the test) checks the test itself.
%
% Results (Greenland nx150 and nx200, Weertman sliding with m=3 and m=1, n=3): order 2.00 (1.99 to 2.04) for all directions and states
% (converged state, velocities x0.5 and x2, 30% smooth perturbation, state after one time step); smallest relative error in the Newton
% direction between 4e-8 and 1e-4; K exactly symmetric. The uv Jacobian is consistent. Check of the test: a term of relative size
% 1e-3 added to the diagonal of K gives order 1.00 and a relative error of 3.8e-3 (flagged), one of 1e-6 still order 1.01 with a
% relative error of 3.8e-6.
%
% Found with the test: CtrlVar.NRviscosity and CtrlVar.NRbeta2 had no effect, as uvMatrixAssemblySSTREAM.m always forms the full
% Newton Jacobian. Removed (17 Oct 2026), together with the also unused CtrlVar.Picard and CtrlVar.Picarditmax, and the parfeval variant
% of the uv assembly (uvMatrixAssemblySSTREAM_Parfeval.m, uvMatrixAssemblySSTREAMpartitionTriplets.m, and its branch and the related
% test block in uvMatrixAssembly.m), the only place where CtrlVar.NRviscosity was still read.
%
% * Mesh key: robust detection of changes of the mesh in UpdateMUA (17 October 2026):
%
% UpdateMUA.m only recognised a changed mesh by changed numbers of nodes, elements and nodes per element. A change of the connectivity
% with the same sizes (eg a renumbering of nodes, or a reordering of elements) or a change of the coordinates went unnoticed, and the
% quantities derived from the mesh were left unchanged. Test with the Greenland nx150 mesh: after a reordering of the elements (same
% sizes and same sum of the connectivity), the uvh Jacobian assembled with the outdated quantities had a relative error of 1.0.
%
%   New ArrayChecksum.m: exact, platform-independent checksum of the bit patterns of a numerical array (two weighted sums modulo primes
%   below 2^31, all intermediate values below 2^53; the probability that a change goes undetected is about 1 in 4e18).
%   New MeshKey.m: Key.Topology (checksum of the connectivity, its size, and the number of nodes) and Key.Geometry (checksum of the
%   coordinates and their size). About 25 ms for the nx150 mesh, about 120 ms for a mesh with 1e6 elements.
%   CreateMUA.m and UpdateMUA.m store the key in MUA.MeshKey. UpdateMUA.m compares it with the key of the current mesh (in addition to the
%   sizes): if the topology or the geometry has changed, the derived quantities are recalculated; the assembly patterns only if the
%   topology has changed. If MUA has no key (eg from an older file), everything is recalculated once.
%   AssemblyPatternCache.m stores the topology key in the patterns (field TopologyKey).
%   uvhMatrixAssemblySSTREAM_SPMD.m: the mesh part of the fingerprint used to decide whether to resend the data to the workers is now
%   the mesh key, instead of [Nnodes Nele sum(connectivity(:))], which did not detect a renumbering or a change of the coordinates.
%   UaPardisoSolve.m: identifies structural patterns by TopologyKey (if available) instead of comparing i0 and j0.
%
% Two errors in UpdateMUA.m found with the tests (both only relevant after a change of the mesh):
%   - The Cholesky factor of the mass matrix (MUA.MC, MUA.Mp; CtrlVar.MUA.CholeskyMassMatrix=true) was not recalculated together with
%     the mass matrix, and was left for the old mesh (relative error 1.8e-3 after moving one node by 50 m). Now recalculated.
%   - The workers (BuildMuaWorkers.m) were built before the element areas were recalculated, so they got the old areas, which are used
%     in the SUPG terms of the uvh assembly (parallel and serial Jacobians differed by 1.3e-5 after moving one node). The element areas
%     are now recalculated before the workers are built.
%
% Tests (nx150): changes detected: local node swap, element swap, renumbering of two nodes, change of one coordinate in the last bit, swap
% across a block boundary. UpdateMUA: unchanged mesh 0.04 s, nothing recalculated; one node moved: Deriv, DetJ and areas changed in
% exactly the 8 elements around it, patterns kept; elements reordered: patterns rebuilt, Jacobian identical to the original (1.5e-16);
% parallel assembly (4 thread workers) identical to the serial one after both changes, data resent only when the mesh changed.
% Transient run, 5 time steps: results bitwise identical to the run before these changes, 86.6 s instead of 88.5 s.
%
% * PARDISO: structural sparsity pattern of the mesh as the fixed pattern (17 October 2026):
%
% Previously UaPardisoSolve.m learned the fixed pattern from the matrices it was given, and extended it whenever a matrix had a new
% nonzero. As sparse() drops entries that are exactly zero, and whether a Jacobian entry is exactly zero depends on the state (ice free,
% floating, current velocities, ...), new nonzeros appeared during the first iterations, and each required a new analysis (eg 10
% extensions in 24 solves in a Greenland nx200 run from scratch, where PARDISO then gained only 5%). Now the structural pattern of the
% mesh, ie all entries the assembled matrices can have, is used: the assembly patterns already stored in MUA.uvAssemblyPattern and
% MUA.uvhAssemblyPattern (computed by AssemblyPatternCache.m in CreateMUA/UpdateMUA, fields i0, j0, neq) give it directly. No changes
% to MUA were needed.
%
%   solveKApe.m, solveKApeSymmetric.m: new optional ninth input StructuralPattern, passed on to ABfgPreEliminate.m (new optional 7th
%   input), KKT_null_space_decomposition_chol.m (new optional 6th input; for the direct reduction for nodal constraints, the pattern of
%   H(f,f) is passed as a structure with fields Pattern and Index) and UaPardisoSolve.m (new optional 6th input).
%   SSTREAM_TransientImplicit.m passes MUA.uvhAssemblyPattern, SSTREAM2dNR2.m passes MUA.uvAssemblyPattern. Other callers are unchanged
%   and use the learned pattern as before.
%   UaPardisoSolve.m: builds the logical pattern (tril for the symmetric types, plus the diagonal) once per mesh and caches it; a state
%   gets it when it is created, or when a different structural pattern arrives (eg after remeshing). The pattern is ignored if it does
%   not match the size of the matrix. If a matrix has an entry outside the pattern, the pattern is extended as before, so that a wrong
%   or outdated pattern only costs time. New statistic StructuralPatterns.
%   New option CtrlVar.LinearSolver.PARDISO.StructuralPattern=true (false: learned pattern, as before).
%
% Tests (Greenland, nx150): with a Jacobian with 2% of its off-diagonal entries zero followed by the full one, 1 analysis instead of 2 for
% each of uvh with constraints, uv without constraints and uv with nodal constraints; solutions agree with MATLAB to 1e-13. A pattern of
% the wrong size is ignored, one that is too small is extended. Transient run, 5 time steps from a restart file, 31 Newton-Raphson
% iterations, 4 thread workers: 1 analysis for 31 factorisations (previously 7), 88.5 s instead of 93.0 s (MATLAB solver: 159.5 s); uvh
% solves 54.9 s instead of 61.4 s; results agree with the MATLAB solver to 1e-12.
%
% * PARDISO (Intel oneMKL) as an alternative linear solver (16 October 2026):
%
% New option CtrlVar.LinearSolver.Type="PARDISO" (default "MATLAB"): the linear systems solved in ABfgPreEliminate.m (uv and uvh
% Newton-Raphson systems) are solved with the multithreaded sparse direct solver PARDISO of Intel oneMKL instead of decomposition/mldivide
% (UMFPACK). In solveKApe.m this covers AsymmSolver="Bempty" (no constraints, solved in solveKApe.m itself) and
% "EliminateBCsSolveSystemDirectly" (nodal constraints, B*B' diagonal, solved in ABfgPreEliminate.m), ie the cases selected with
% AsymmSolver="auto"; with PARDISO, no decomposition is returned (dA=[]), and a given dA (modified Newton-Raphson) gives an error.
% Bempty test (Greenland, 168,210 unknowns): 0.25 s instead of 2.6 s per solve.
%
% Symmetric systems (solveKApeSymmetric.m): with PARDISO, SymmSolver="Bempty" is solved as symmetric indefinite (PARDISO matrix type -2,
% LDL'), and in "NullSpace" (KKT_null_space_decomposition_chol.m, new optional fifth input CtrlVar) the reduced Hessian Hr=Z'*H*Z is
% factorised with a Cholesky factorisation (type 2); if that fails (Hr not positive definite) the system is solved as symmetric
% indefinite instead, with a warning. For the symmetric types only tril(A) is passed. Tests with the uv system of the Greenland model
% (112,140 unknowns): Bempty 0.12 s instead of 0.51 s; NullSpace with 2324 constraints 0.27 s instead of 0.61 s (the Cholesky
% factorisation of Hr 0.04 s instead of 0.41 s; the QR factorisation of A' and Z'*H*Z then took most of the time).
%
% KKT_null_space_decomposition_chol.m: new direct route for nodal constraints (each row of A has one nonzero, d_k in column c_k, and
% all columns are different). Then Z=I(:,f), with f the free unknowns, and x(c)=b./d, H(f,f) x(f) = g(f)-H(f,c) x(c),
% lambda=(g(c)-H(c,:) x)./d. Identical to the general null-space solution, but without the QR factorisation of A', Z'*H*Z, and the
% factorisation of A*A'. For all other constraints (eg ties between unknowns, even if A*A' is diagonal) the general route is used.
% uv system (112,140 unknowns, 2324 nodal constraints): 0.49 s instead of 0.61 s with MATLAB, 0.16 s instead of 0.27 s with PARDISO;
% results identical to the general route.
%
% UaPardisoSolve.m: new optional fifth input MatrixType ("nonsymmetric", "spd", "symmetric"); separate states for each size and type,
% at most 6 (the least recently used one is released); systems with fewer than CtrlVar.LinearSolver.PARDISO.MinimumSize (default 20000)
% unknowns are solved with mldivide, and diagonal systems directly (for these, the analysis by PARDISO would cost more than it saves).
%
% solveKApeSymmetric.m: corrected isfield(CtrlVar,"Symmsolver") to isfield(CtrlVar,"SymmSolver"). Previously CtrlVar.SymmSolver was
% always set to "auto", so that a solver given in CtrlVar.SymmSolver was ignored.
%
% Linux (NOT YET TESTED WITHIN MATLAB): pardisoMex.c now also compiles on Linux (dlopen instead of LoadLibrary; BuildPardisoMex adds -ldl).
% The library is libmkl_rt.so.N (N=3 for MKL 2026, N=2 for MKL 2021-2025), eg .../lib/libmkl_rt.so.3 after pip install mkl. MKL is loaded
% with RTLD_LOCAL|RTLD_DEEPBIND, so that it does not use symbols of the MKL loaded by MATLAB, and libmkl_core, libmkl_intel_lp64/ilp64 and
% the threading layer are loaded first, from the folder of libmkl_rt (they have no RUNPATH, and LD_LIBRARY_PATH can not be changed within
% MATLAB). Tested outside MATLAB (Ubuntu 24.04, gcc 13, MKL 2026.2, with a stand-in for the MEX API): all commands, both threading layers.
% Not tested: the interaction with the MKL, TBB and OpenMP libraries of MATLAB on Linux. macOS is not supported.
%
% New folder UaSource/UaPardiso with:
%   pardisoMex.c, BuildPardisoMex.m  MEX interface to PARDISO (compile with BuildPardisoMex, eg with MinGW64)
%   UaPardisoSolve.m                 manages the PARDISO states, the fixed sparsity patterns, a residual check, and a fallback to decomposition
% New settings: CtrlVar.LinearSolver.PARDISO.Library (path to mkl_rt, eg ...\Library\bin\mkl_rt.3.dll; if empty the environment variable
% UA_MKL_RT is used), .ThreadingLayer ("TBB" or "SEQUENTIAL") and .ResidualTolerance (1e-10).
%
% Notes on the implementation:
%  - The MKL runtime library is not part of MATLAB (the MKL in MATLAB does not contain PARDISO) and must be installed separately, eg
%    pip install mkl. It is loaded at run time by pardisoMex, so no MKL headers or import libraries are needed to compile it.
%  - MKL's OpenMP threading layer can not be used within MATLAB: it requires a newer OpenMP runtime than the one MATLAB has loaded, and
%    MKL then aborts the MATLAB process. pardisoMex therefore sets the TBB threading layer (or the sequential one) itself, before the
%    first computation, and rejects the OpenMP layer. With the TBB layer MKL uses all cores.
%  - The matrix is passed without copying: the CSC arrays of a MATLAB sparse matrix A (64-bit indices) are the CSR arrays of A', which
%    are passed to the 64-bit interface pardiso_64 with zero-based indexing, and PARDISO solves the transposed system (iparm(12)=2).
%  - Fixed pattern: the analysis (reordering, symbolic factorisation) is much more expensive than the numerical factorisation, but the
%    pattern of the Newton-Raphson matrices changes from one iteration to the next (exact zeros are not stored in MATLAB sparse matrices,
%    and the rows of constrained nodes change with the active set). UaPardisoSolve therefore keeps, for each matrix size, a pattern that
%    contains all patterns seen so far (initially the pattern of the full Jacobian and the diagonal); pardisoMex copies the values into
%    this pattern, with explicit zeros, and the analysis is only repeated if the pattern has to be extended.
%
% Results (Greenland, 56,070 nodes, Atilde with 168,210 unknowns, 16 cores): factorisation 0.135-0.18 s with a stored analysis and
% solve 0.09 s, against 2.5 s and 0.17 s with decomposition (UMFPACK, which does not benefit from more threads for these matrices);
% analysis 1.2-1.4 s. solveKApe 0.40 s instead of 2.6 s. Model run of 5 time steps with the parallel assembly (4 thread workers): 93.0 s
% instead of 159.5 s with decomposition, with the same Newton-Raphson iterations and results agreeing to 1e-12 (7 analyses for 31
% factorisations); the uvh solves took 61 s instead of 128 s.
%
% * New function uvhAssemblyScalingTest.m: performance of the parallel uvh assembly on a given computer (16 October 2026):
%
% For a model state, read from a file saved with CtrlVar.uvhSaveNewtonSystem.Save=true or from a restart file, the function times the
% serial uvh assembly with different numbers of threads, the solution of the linear system, and the parallel (spmd) assembly with
% different numbers of workers (default 1, 2, 4, 8, 12, ... up to the number of cores). It fits the model T(n)=a+b/n+d*n to the times of
% the parallel assembly, and recommends either the parallel assembly with the smallest number of workers for which the time is within 5%
% of the smallest time, or the serial assembly if the parallel one is not faster. It also estimates the reduction of the time per
% Newton-Raphson iteration. Parallel pools opened by the function are closed at the end, and a pool that was open before is restored.
%   Results=uvhAssemblyScalingTest("RestartFile.mat") ;
%
% Tests: Greenland (25,095 nodes, 16 cores): 8 workers recommended, assembly 2.7 times faster than the best serial assembly, estimated
% reduction of the time per Newton-Raphson iteration 27%. Thule (5,995 nodes): serial assembly recommended (the parallel assembly of the
% residual only is slower than the serial one, and the full assembly only slightly faster).
%
% Larger Greenland mesh (56,070 nodes, 110,976 elements, 16 cores): 9 or 10 workers recommended in repeated runs, assembly 2.35 times
% faster than the best serial assembly (1.52 s, 8 threads), estimated reduction of the time per Newton-Raphson iteration 22%. Model
% parameters a=0.231 s, b=2.507 s, d=0.0152 s, optimum at about 13 workers; the time is within 7% of its minimum for 8 to 16 workers.
%
% uvhAssemblyScalingTest.m now also supports thread-based pools: Options.Profile="Threads" (parcluster is not supported for
% thread-based pools, and the profile of an open pool is now determined with a local function PoolProfile). Results for the larger
% Greenland mesh with a thread-based pool: the transfer of the Jacobian to the client is about 10 times cheaper than with processes
% (0.01 s instead of 0.13-0.17 s), and with few workers the thread workers are faster (they use implicit multithreading), but beyond
% about 8 workers they compete for the cores of the one process and the time increases again (1 worker 1.81 s, 4: 0.75 s, 8: 0.70 s,
% 16: 1.29 s; model a=0.072 s, b=1.732 s, d=0.0658 s). 4 workers recommended, estimated reduction of the time per Newton-Raphson
% iteration 21% (process-based pool: 9 or 10 workers, 22%).
%
% Also checked: MATLAB decomposition objects can not be passed to or from the workers of either process-based or thread-based pools
% (they are serialised, and arrive as the decomposition of an empty matrix), but they can be computed and used on the workers, eg within
% an spmd block (kept as a Composite), or with parallel.pool.Constant(@() decomposition(A)). The parallel uvh assembly is not affected.
% See the section "Process-based and thread-based pools" of the chapter on the parallel assembly in the Ua Compendium.
%
% uvhMatrixAssemblySSTREAM_SPMD.m has a new optional fifth output, Timing, with the times of the parts of the parallel assembly
% (fingerprint, sending changed inputs, the spmd block, the longest time on a worker, the transfer of the results, total). It is used by
% uvhAssemblyScalingTest.m. Without the fifth output the behaviour is unchanged.
%
% * Parallel (spmd) uvh assembly: overheads reduced (16 October 2026):
%
% uvhMatrixAssemblySSTREAM_SPMD.m (used if CtrlVar.Parallel.uvhAssembly.spmd.isOn=true and a parallel pool is open) has been revised.
% Measurements (Greenland, 25095 nodes, 8 workers) had shown that the largest costs of the previous version were not the assembly on the
% workers (about 0.19 s), but sending the inputs to the workers in each call (about 0.21 s, mainly F0 and F1 complete, and CtrlVar, a
% large nested structure of only 0.3 MB that alone took about 0.13 s), so that the parallel assembly was hardly faster than the serial
% one (0.71 against 0.78 s), and residual-only assemblies were slower (0.36 against 0.23 s). Changes:
%
% (1) F0, F1, CtrlVar (except CtrlVar.uvhMatrixAssembly), UserVar and BCs1 are kept on the workers between calls. Only ub, vb, h, b and s,
% which change within a uvh solve, and CtrlVar.uvhMatrixAssembly are sent in each call. Everything else is only sent again if it has
% changed: F0 and F1 are compared with a fingerprint (weighted sums of all numerical fields except the changing ones, the mesh, the number
% of workers and the pool), CtrlVar, UserVar and BCs1 with isequaln. Whatever changes, the worst case is an unnecessary resend.
% (2) The assembly and the summation over the workers are done in one spmd block instead of two.
% (3) BuildMuaWorkers.m now adds the sparsity pattern of each sub-domain to the worker data, so that the workers use the faster
% pattern-based construction of the sparse matrix, as the serial assembly does.
%
% Results (25095 nodes, 8 workers, with another run using 6 cores in the background): full assembly 0.38 s (serial 0.78 s, previous
% parallel version 0.71 s), residual-only assembly 0.19 s (serial 0.23 s, previous 0.36 s), complete resend about 0.27 s. In a run of 10
% time steps the wall time was reduced by 18% compared with the serial assembly, with the same Newton-Raphson iterations and results
% agreeing to 1e-12 (round-off, different summation order). Tested also with remeshing (Thule, 5 remeshings in 20 time steps): identical
% meshes and results agreeing to round-off. For small meshes the parallel assembly does not pay off: Thule with about 12,000 elements and
% 4 workers was 10% slower than with the serial assembly.
%
% Scaling with the number of workers (otherwise idle computer with 16 cores, same Greenland model): full assembly 1.25 s with 1 worker,
% 0.64 s (2), 0.44 s (4), 0.32 s (8), 0.31 s (12) and 0.35 s (16); residual only 0.48 s (1) to 0.17 s (8). The serial assembly benefits
% from implicit multithreading only up to about 4 threads (1.12 s with 1 thread, 0.73 s with 4, no further gain with 8 or 16). The
% measured times are described by T(n)=a+b/n+d*n with a=0.073 s (mainly the transfer of the Jacobian to the client), b=1.154 s (the
% assembly, divided between the workers) and d=0.0126 s (data sent to each worker, summation over the workers, competition for cores
% and memory bandwidth), with an optimum at about 10 workers (speed-up of about 2.3 compared with the best serial assembly). See the
% chapter "Parallel assembly of the uvh system: performance and scaling" in the Ua Compendium.
%
% Possible further improvements (not implemented): sending each worker only the values at the nodes of its own elements, so that the
% total amount of data sent no longer grows with the number of workers, and collecting the contributions of the workers on a fixed
% global sparsity pattern instead of summing sparse matrices. Estimated to reduce the full assembly with the optimal number of workers
% from about 0.31 s to about 0.2-0.25 s, and to shift the optimum to more workers; about 5% of the time of a Newton-Raphson iteration on
% a computer with 16 cores.
%
% * Testing the Jacobian of the uvh system (15 October 2026):
%
% (1) New function TestuvhJacobian.m: a Taylor test of the uvh Jacobian. For a state x with residual R and Jacobian K (from uvhAssembly.m)
% and a direction d, the remainder ||R(x+epsilon d)-R(x)-epsilon K d|| is calculated for epsilon=1e-1,...,1e-9, separately for the u, v
% and h equations. For a consistent Jacobian it decreases as epsilon^2 (order 2), for a missing or wrong term as epsilon (order 1). Directions:
% the Newton direction of the iteration, and smooth random perturbations of u, v, h, and of all three. Two measures are reported: the
% observed order above the round-off level, and the smallest relative error over all epsilon (a missing term of relative size delta gives
% at least delta). A block is flagged if the order is below 1.5 and the smallest relative error above 1e-4. Usage:
%   Results=TestuvhJacobian(FileName)                                 % all states saved in a file (see (2))
%   Results=TestuvhJacobian(FileName,Options)                         % eg Options.States, .Directions, .CtrlVar, .Plot
%   Results=TestuvhJacobian(UserVar,CtrlVar,MUA,F0,F1,l1,BCs1,Options)  % a single state
%
% (2) New functions SaveuvhNewtonState.m and uvhNewtonStateSelected.m, and new options: CtrlVar.uvhSaveNewtonSystem.Save=true saves the
% state of selected Newton-Raphson iterations of the uvh solver (CtrlVar.uvhSaveNewtonSystem.Iterations, .RunSteps, .MaxStates) to one
% file (CtrlVar.uvhSaveNewtonSystem.FileName), each state as a variable State_0001, State_0002, ..., with the mesh saved only once (and
% again after remeshing). CtrlVar.uvhTestJacobian.Test=true runs TestuvhJacobian.m directly during a run, in selected iterations. Both
% are off by default (no change in behaviour). Called in SSTREAM_TransientImplicit.m after the Newton direction has been calculated.
%
% (3) Result: the Jacobian is consistent. RadialIceCap (3 states), Greenland with 6324 nodes (21 states, 3 time steps, Newton-Raphson
% iterations 1-5) and Thule with a calving front and grounding lines (24 states): for all random directions and all blocks the observed
% order was 1.98-2.03. For the Newton direction, the order was below 1.5 in some blocks in late iterations, but only where its
% second-order remainder is smaller than the round-off noise of the residual (mainly the h equation for thick ice, where the terms are
% large); the smallest relative errors there were 1e-6 to 1e-10. No block was flagged. The phase of slower-than-quadratic convergence
% of the u and v residuals in the uvh solve (reductions by a factor of about 10 per iteration with full steps, after the h residual has
% been removed in the first iteration) is therefore not caused by an inconsistent Jacobian, but presumably by the nonlinearity of the
% flow and sliding laws (pre-asymptotic phase of the Newton-Raphson method).
%
% * Reuse of the Jacobian (modified Newton) in the uvh solver: implemented, tested, and found not to pay off (15 October 2026):
%
% (1) solveKApe.m can now return and reuse the matrix decomposition: [x,y,dA]=solveKApe(A,B,f,g,x0,y0,CtrlVar) returns the decomposition,
% and [x,y]=solveKApe(...,CtrlVar,dA) uses it without a new factorisation. Supported for "Bempty" (no constraints; decomposition() instead
% of backslash) and "EliminateBCsSolveSystemDirectly" (point constraints, via ABfgPreEliminate.m, which no longer builds the eliminated
% matrix if the decomposition is given). Calls with two outputs and without dA are unchanged (backslash for "Bempty").
%
% (2) SSTREAM_TransientImplicit.m: optional reuse of the Jacobian and its decomposition within a uvh solve. The Jacobian is reused if the
% previous step was a full Newton step that reduced the residual by at least the factor CtrlVar.ModifiedNRuvhReductionCriterion (now
% default 1e-2), at most CtrlVar.ModifiedNRuvhIntervalCriterion-1 times in a row. A step with the reused Jacobian is only tried at full
% length, and accepted if it reduces the residual by at least the factor CtrlVar.ModifiedNRuvhAcceptCriterion (new, default 0.1) or meets
% the convergence criteria; otherwise it is discarded and a full Newton step is done at the same point. Default
% CtrlVar.ModifiedNRuvhIntervalCriterion=1, ie off, which gives exactly the previous behaviour (tested: bit-identical Newton-Raphson paths
% and results).
%
% (3) Result: reusing the Jacobian did not pay off (Greenland, 25095 nodes, 10 time steps). Reuse in the final iterations (reduction
% criterion 1e-2): 73 instead of 61 Newton-Raphson iterations, 60 instead of 61 factorisations, run time +6% (interval 2) and +10%
% (interval 3). Reuse after every full step (reduction criterion 1, acceptance criterion 0.5): 101 and 108 iterations, 55 and 75
% factorisations, run time +8% and +34%. The steps with the reused Jacobian were accepted, but reduced the residual much less than full
% Newton steps (typically by a factor of about 10, against full steps that reach the convergence tolerance directly), and so added
% iterations rather than replacing full ones. This agrees with earlier experience with the classical modified Newton method.
%
% * Line search of the uv and uvh Newton-Raphson solvers (15 October 2026):
%
% (1) New line search rLineminUa3.m, now the default (CtrlVar.NewtonLineSearch.Method="rLineminUa3"). The logic is identical to that of
% the previous rLineminUa2.m, but unnecessary evaluations of the residual (cost function) have been removed: (a) In
% SSTREAM_TransientImplicit.m the residual at gamma=0 was evaluated twice, although it had just been assembled together with the
% Jacobian. It is now calculated from the assembled residual (uvhResidualMeasures.m). (b) The residual at the full Newton step was
% evaluated twice, once before the line search, and again at the end of rLineminUa2.m to obtain rForce, rWork, D2 and rBlocks.
% rLineminUa3.m keeps all outputs of each evaluation, and evaluates the cost function as func(1,gamma*s) rather than func(gamma,s), which
% is arithmetically identical. The Newton-Raphson paths and results are bit-identical to those with rLineminUa2.m. Tests: 70-72% fewer
% residual evaluations in the uvh solver (Greenland with 6324 and 25095 nodes, Thule), and 18% shorter run times (Greenland, 25095 nodes,
% 10 time steps, about 0.6 s less per Newton-Raphson iteration). In a 10-year run (Greenland, 25095 nodes, automated time stepping,
% 97 time steps, 654 Newton-Raphson iterations), the run time was reduced by 25%, from 1851 s to 1390 s, with 810 instead of 2749 residual
% evaluations (0.71 s less per Newton-Raphson iteration), and with bit-identical results. The uv solver (SSTREAM2dNR2.m) also uses
% rLineminUa3.m by default, with results identical to those of the previous combination of a full-step test and rLineminUa.m.
% Correction (15 October 2026, after the first version of rLineminUa3.m had been committed): in the Cauchy step of the fall-back, the
% step length can be a sparse scalar, which caused an error in the look-up of stored evaluations (sprintf does not accept sparse
% input). The step length is now converted with full() for the look-up only. Found in the Thule example, in the initial uv solve.
%
% (2) CtrlVar.NewtonLineSearch.Method selects the line search of the uv and uvh solvers: "rLineminUa3" (default), "NewtonLineSearch"
% (see (3)), or "rLineminUa2" (the previous implementation, including the unnecessary evaluations; for comparison only).
%
% (3) New function NewtonLineSearch.m (optional, CtrlVar.NewtonLineSearch.Method="NewtonLineSearch"): backtracking with the Armijo
% condition, quadratic and cubic interpolation with safeguards [0.1,0.5], an optional refinement after backtracking (minimum of a
% parabola, repeated while it halves the residual), and optional extrapolation. Parameters: CtrlVar.NewtonLineSearch.ArmijoParameter
% (default 0.25), GammaMin (1e-3), MaxBacktrackingSteps (15), RefinementSteps (2), RefinementRatio (0.5), Extrapolation (false),
% MaxExtrapolationSteps (3) and ExtrapolationFactor (2). In tests it required fewer residual evaluations than rLineminUa3.m, but more
% Newton-Raphson iterations (Greenland, 6324 nodes, 30 time steps: 176 against 150; Thule, 20 time steps: 232 against 229; initial uv
% solve: 15 against 12). With the textbook value ArmijoParameter=1e-4 the number of iterations was considerably larger (193, 266 and 23),
% as full steps that reduce the residual only marginally are then accepted, although a shorter step often reduces it much more. As a
% residual evaluation is cheap compared with a Newton-Raphson iteration (about 1:12 in these tests), additional evaluations along the
% Newton direction that find a better step are worthwhile, which is what the logic of rLineminUa2/3.m does. NewtonLineSearch.m is
% therefore not the default.
%
% (4) The option to use the work residual as the merit function of the uv and uvh solvers, CtrlVar.uvhMinimisationQuantity and
% CtrlVar.uvMinimisationQuantity="Work Residuals", has been removed. The force residual is always used. The work residual is still
% calculated, and used in the convergence criteria. (Changes in CalcCostFunctionNRuvh.m, CalcCostFunctionNR.m, JuvhLSQ.m and
% Ua2D_DefaultParameters.m. CtrlVarValidityCheck.m prints a note, and removes the fields, if these options are still set.)
%
% (5) The calculation of the residual measures (rForce, rWork, D2 and rBlocks) has been moved, unchanged, from CalcCostFunctionNRuvh.m
% and CalcCostFunctionNR.m to the new functions uvhResidualMeasures.m and uvResidualMeasures.m, so that these measures can also be
% calculated from an already assembled residual. (In CalcCostFunctionNRuvh.m, rBlocks is now always calculated.)
%
% (6) Bug fix in BackTracking.m: if additional outputs of the function were requested (nOut>1), these were from the last function
% evaluation, which is not necessarily the returned (best) step. They are now from the returned step, with a further evaluation if
% needed. In tests with synthetic functions this affected about 6% of the cases. This changes the results of callers that use these
% outputs, eg MassContinuityEquationNewtonRaphson.m, wherever the returned step was not the last one evaluated.
%
% * Grounding-line plots now based on h-h_{gl}, ie thickness minus grounding line thickness, rather than the flotation mask.
% Thanks to Deborah Rhee for this suggestion. The grounding lines may now look a bit smoother. (Note that the grounding line
% position itself is never used in the code when finding the solutions.)
%
% * Default for the positive-thickness constraints: active set on, and thickness penalty off.
%
%   CtrlVar.ThicknessConstraints=true;
%   CtrlVar.ThicknessPenalty=false;
%
% To use the penalty in combination with the active set, set
%
%   CtrlVar.ThicknessConstraints=true;
%   CtrlVar.ThicknessPenalty=true;
%
% * An error in the calculation of the adjoint gradient for Budd law was spotted and corrected. In one line of the code, a term was
% multiplied twice by the grounding/floating mask, instead of once. This will not have caused any errors in the forward
% solution, but might have slowed done convergence in an inversion.
%
% * Thickness penalty (SoftPlus option), changes made 8-9 October 2026:
%
% (1) K is again a mass-balance rate (units 1/time). Between 31 August and 9 October 2026 K was internally divided by the time step,
% making the penalty a per-time-step term. This was reverted because the penalty then changed the thickness by a dt-independent
% amount in every time step, so that the solution had no limit for dt->0, the number of Newton iterations did not decrease for
% small dt, and the automated time stepping could drive dt to very small values. With a rate, a node is held above ThickMin by the
% penalty alone if theta*K*delta > |a|max, otherwise the active set takes over.
%
% (2) The penalty is centred at hmin=ThickMin+delta with delta=max(ThickMin,deltaAbs), default deltaAbs=0.1. This is the previous
% hmin=2*ThickMin for ThickMin>=deltaAbs, but does not collapse to zero for ThickMin=0.
%
% (3) By default the smoothing distance l is calculated at runtime as l=lRel*delta (default lRel=0.1), ie from the value of ThickMin
% in the input file. An explicitly set l is used as given.
%
% (4) The penalty can be applied as a lumped nodal term, CtrlVar.ThicknessPenaltyNodal=true (default false), which excludes nodes with
% thickness boundary conditions or active thickness constraints.
%
%   CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.K=100;     % 1/time
%   CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.l=NaN;     % l=lRel*delta
%   CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.lRel=0.1;
%   CtrlVar.ThicknessPenaltyMassBalanceFeedbackSoftPlus.deltaAbs=0.1;
%
% * Active-set method for the positive-thickness constraints, changes made 8-9 October 2026:
%
% (1) Constraints are released based on the sign of the discrete reaction at the node (CtrlVar.ActiveSet.ReleaseCriterion="discrete").
% The previous approach, based on reactions mapped to a nodal field with the inverse of the consistent mass matrix, introduced sign
% oscillations at individual nodes and caused the active set to cycle. It can still be selected with ReleaseCriterion="consistent-mass".
%
% (2) If a node released in an active-set iteration is again below ThickMin after the following solve, it is kept constrained, and no
% further nodes are released for the rest of that time step. This prevents cycling. Re-activation of such nodes is not counted as a cycle.
%
% (3) If the active set is nevertheless found to be cyclical, one final uvh solve is done with the last active set, so that the returned
% solution is consistent with the returned active set.
%
% (4) Initial active set: nodes with h<=ThickMin (rather than h<ThickMin) are included, so that a previous active set is recovered after
% restarts and re-meshing.
%
% (5) CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints has been removed. It caused small updates of the active set to be
% discarded, which could leave thickness violations in place over several time steps.
%
% (6) The active set is returned from the semi-implicit uv-h solver and carried over to the next time step.
%
% * Explicit estimate at the start of each time step (CtrlVar.ExplicitEstimationMethod="-Adams-Bashforth-"), corrected 9 October 2026:
%
% The rates used in the explicit estimate (F.dhdt, F.dubdt, F.dvbdt, F.duddt, F.dvddt) are backward differences over the previous time
% step (see UpdateFtimeDerivatives.m), ie the time derivatives at the midpoint of that step. These were previously used in the
% variable time step Adams-Bashforth formula (ExplicitEstimation.m), which requires the time derivatives at the start of the steps.
% As a result the estimate was only first-order accurate (for a constant time step it reduced to X+dt*dX/dt). The estimate is now
% calculated with ExplicitEstimationUsingBackwardDifferences.m, which is second-order accurate also for variable time steps, and exact
% if the fields vary quadratically in time. The name "-Adams-Bashforth-" is kept for compatibility. The estimated thickness is not
% allowed to be below ThickMin. Related changes:
%
%   - UaFields has a new property dtRates, the time step over which the rates were calculated.
%   - Fm1 (the rates of the previous time step) is a UaFields object, is written to restart files, and is mapped onto a new mesh
%     after remeshing (MapFm1BetweenMeshes.m). For older restart files Fm1 is created with all rates set to NaN, and the estimate then
%     falls back to linear, or no, extrapolation in the first time step(s).
%   - Which level of extrapolation (second order, linear, none) can be used is determined from the available data, node by node, and
%     no longer from the run-step number.
%   - The time step is taken from F.dt (and not from CtrlVar.dt) in UpdateFtimeDerivatives.m and ExplicitEstimationForUaFields.m.
%     ExplicitEstimationForUaFields.m has the new time step as an additional (11th) input argument.
%
%
% * Reactions of the thickness constraints as an output (15 October 2026):
%
% The equivalent local mass-balance rate of the thickness constraints is now provided as an output field, F.ah (units of distance per
% time, eg m/yr), positive where the constraints add ice and negative where they remove ice. It is calculated at the end of each time
% step from the discrete reactions of the constraints, and includes all thickness constraints: the active set (h>=ThickMin), the
% level-set constraints downstream of calving fronts, and user-defined thickness boundary conditions. In addition, the total rates at
% which mass is added and removed by the constraints (units of mass per time, eg kg/yr) are stored for each time step in
% RunInfo.Forward.ThicknessConstraintsMassAdded and RunInfo.Forward.ThicknessConstraintsMassRemoved. Downstream of a calving front, for
% example, the mass removed is the calving flux. The mass budget of each time step closes exactly, ie the change in ice mass is equal to
% the integrated mass balance, minus the net outflow across the boundary of the domain, plus the mass added and minus the mass removed by
% the constraints (tested to round-off in the Thule example). The constraints can make a substantial contribution to mass budgets: in a
% model of the Greenland ice sheet they added about 100 Gt/yr at ice-free nodes along the margins, as the thickness equation of such a
% node also contains a share of the mass balance and thinning of neighbouring ice-covered nodes (through the overlap of the shape
% functions). The calculation is done in the new function ThicknessConstraintReactions.m, which is also used by ActiveSetUpdate.m.
%
% * Time-discretisation error estimate near calving fronts (15 October 2026):
%
% The time-discretisation error estimate (TimeDiscretisationErrorEstimate.m) is now also masked at all nodes of elements cut by the
% calving front (F.LSFMask.NodesOn), and not only at the nodes downstream of the front (F.LSFMask.NodesOut). The thickness at the nodes
% of cut elements depends on the position of the front within the element, and does not evolve smoothly in time.
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
% active-set iterations. The change in results is at most ReleaseAlpha*ThickMin (ie 1 mm for ThickMin=1 m) per time step at such nodes.
%
% Each of the two changes on its own did not reliably reduce run times, but together they reduce the number of uvh solves by about a
% third, as the active set is then more often correct after the first solve. Tests (Greenland, 6324 and 25095 nodes, 10 years, automated
% time stepping): 16-19% fewer Newton-Raphson iterations and 17-18% shorter run times, with no loss of accuracy compared with
% high-accuracy reference solutions. For a synthetic ice cap (RadialIceCap) and for the Thule example with a moving calving front and
% mesh refinement, the changes had practically no effect on run times or results. To obtain the previous behaviour, set
% CtrlVar.ActiveSet.PredictiveActivation=false and CtrlVar.ActiveSet.ReleaseAlpha=0. Predictive activation should only be used together
% with ReleaseAlpha>0; a note is printed otherwise.
%
% A third, related option, CtrlVar.ActiveSet.ActivationTolerance (default 1e-3), is the counterpart of ReleaseAlpha for activations: a
% free node is only activated if h<ThickMin*(1-ActivationTolerance), ie nodes slightly below ThickMin are accepted as they are (they are not
% set to ThickMin, and mass is conserved). ReleaseAlpha is thus a tolerance on dual feasibility (the sign of the multipliers), and
% ActivationTolerance a tolerance on primal feasibility (the constraint h>=ThickMin), both expressed as a fraction of ThickMin, and with the
% same default value. In the Greenland tests this reduced the number of Newton-Raphson iterations by a further 1-5%, with no measurable
% effect on the results, and no free nodes remained below ThickMin. To obtain the previous behaviour, also set
% CtrlVar.ActiveSet.ActivationTolerance=0.
%
% Related change: The explicit estimate of the thickness at the end of the time step (ExplicitEstimationForUaFields.m) is no longer
% clipped at ThickMin. The (unclipped) estimate is used for predictive activation, and the thickness is made feasible (h>=ThickMin) at the
% start of the uvh solve (uvhRootFinding.m, uvhSemiImplicit.m). With predictive activation switched off, results are unchanged.
%
% * Active-set method: early exit from the uvh iteration if free nodes have non-positive thickness (14 October 2026):
%
% In transient runs with thickness constraints, the uvh iteration is now ended early if nodes that are not constrained have
% non-positive ice thickness in two consecutive Newton-Raphson iterations. The active set is then updated (these nodes are activated),
% and the iteration continued with the same time step. This typically happens after the active-set update has released a group of
% neighbouring nodes, and the unconstrained solution then has thicknesses well below ThickMin at some of them. The momentum equations are
% then locally ill-posed (the viscous stiffness is proportional to h, and the driving stress is not monotone in h), and the
% Newton-Raphson iteration previously often stagnated, and the solution was then accepted at the (looser) acceptable tolerances, or the
% time step reduced. The result of such a solve is discarded anyway, as these nodes are activated in the next active-set update. The
% early exit is only done if a further active-set iteration is possible, so that the result of a time step is always a converged solve.
%
%   CtrlVar.ActiveSet.EarlyExitOnNegativeThickness=true;   % default, set to false for the previous behaviour
%   CtrlVar.ActiveSet.EarlyExitThickness=0;                 % thickness threshold
%
% * Thickness constraints, level set, and related corrections (13-14 October 2026):
%
%   - Duplicate thickness constraints: Nodes with user-defined thickness boundary conditions (BCs.hFixedNode) could also be included in
%     the active set, if they were downstream of a calving front (level-set constraints), resulting in duplicate rows in the constraint
%     matrix. Such nodes are now removed from the active set (ActiveSetInitialisation.m, ActiveSetUpdate.m). In addition, BCs2MLC.m now
%     stops with an error if a node appears more than once in the thickness constraints (BCs.hFixedNode and BCs.hPosNode), also if a
%     node is included more than once in BCs.hFixedNode.
%   - CalculateReactions.m: with the option hReactionsExcludeActiveSet, the wrong Lagrange multipliers were used if there were tied
%     thickness constraints (not used by any current calls).
%   - CtrlVar.LevelSetMinIceThickness has been removed. Downstream of calving fronts the ice thickness is set to CtrlVar.ThickMin, which
%     was already the value used in the active-set method. If CtrlVar.LevelSetMinIceThickness is still set, it is ignored, and a note is
%     printed.
%   - CtrlVar.LevelSetEvolution is now compared case-independently. Previously "-prescribed-" (the default) and "-Prescribed-" were
%     treated differently in some places, so that, for example, a prescribed level set could be mapped instead of being redefined
%     through DefineCalving.m after remeshing.
%   - The information printed at the start of a run (eg "Time-dependent restart run") could be wrong if CtrlVar.TimeDependentRun was
%     set through CtrlVar.ForwardTimeIntegration. This is now based on CtrlVar.TimeDependentRun (CtrlVarValidityCheck.m, PrintRunInfo.m).
%
% * Explicit estimate at the start of each time step: new option, and new default (12 October 2026):
%
% CtrlVar.ExplicitEstimationMethod="-Adams-Bashforth-thickness-only-" is a new option, and is now the default. Only the ice thickness is
% extrapolated (second order, as in "-Adams-Bashforth-"), while the velocities are set to their values at the start of the time step, ie to
% the converged solution of the previous time step. Tests showed that with "-Adams-Bashforth-" the initial residual of the non-linear
% iterations is much smaller than without extrapolation, but that this did not reduce the number of non-linear iterations in realistic
% settings (Greenland): with the same time steps, "-Adams-Bashforth-" required 6-13% MORE iterations than "-no extrapolation-". The reason
% is that the extrapolated velocities are not consistent with the extrapolated thickness (they do not satisfy the momentum equations), and
% velocities can change non-smoothly in time near grounding lines, at fast-flowing outlets, and where the active set of the thickness
% constraints changes. The previous velocities, on the other hand, satisfy the momentum equations for the thickness at the start of the
% time step, and the first Newton-Raphson step from this consistent state is very effective. 
%
%
% * The default value of CtrlVar.ATSTargetIterations has been changed from 4 to 7 (10 October 2026). Runs that do not set this parameter
% will therefore generally use larger time steps than before. In tests (Greenland, nx100, 10 years), ATSTargetIterations=7 gave about 40%
% shorter run times than 4, with differences in ice thickness of at most about 0.5 m (rms below 0.01 m). To obtain the previous behaviour, set
% CtrlVar.ATSTargetIterations=4.
%
% * Restart runs, changes made 10 October 2026:
%
% (1) The active set of the thickness constraints (BCs.hPosNode) and the Lagrange multipliers are restored from the restart file (provided
% the mesh has not changed, and the boundary conditions defined in DefineBoundaryConditions.m are the same as in the previous run). Previously
% the active set had to be recovered in the first time step, and the first uvh solve started with zero multipliers. A restart run now
% continues where the previous run ended, and restarted and uninterrupted runs agree to within round-off errors.
% (2) The run-step number written to the final restart file is now the number of the last completed run step. Previously it was one too
% large, and the restart run skipped one run-step number.
%
% * Estimate of the time-discretisation error of the ice thickness (10 October 2026):
%
% If CtrlVar.TimeDiscretisationErrorEstimate.Use=true, the local time-discretisation error of h is estimated in each time step
% (TimeDiscretisationErrorEstimate.m). This is a predictor-corrector estimate, comparing the second-order explicit estimate of h with the
% computed (implicit, theta=0.5) solution, and is only available for theta=0.5. The estimate is in units of h (ie m), and is the error
% committed in that time step. It is stored as a nodal field in F.hTimeDiscretisationErrorEstimate, and norms of it are stored as time series
% in RunInfo.Forward (hTimeDiscretisationErrorRMS, ..Max, ..Volume, ..ValidArea, etc.). Norms are calculated as sqrt(x'*M*x/A) with the
% consistent mass matrix M, where A is the area over which the estimate is defined. The estimate is not defined (NaN) at nodes with thickness
% boundary conditions, active thickness constraints, h<=ThickMin, and next to nodes whose constraint state changed within the last three time
% steps, as the thickness then does not evolve smoothly in time. Optionally, the estimate can be accumulated over a run
% (CtrlVar.TimeDiscretisationErrorEstimate.Accumulate=true, see F.hTimeDiscretisationErrorAccumulated). The accumulated estimate is an
% "error budget", ie the sum of the local errors, and not the global error, as errors are transported and may decay.
%
% Experience from tests: The estimate is useful as a diagnostic for where and when time-discretisation errors arise. Errors are typically
% largest during initial adjustments of the model state, and at fast-flowing outlets, where the spatial resolution may also be insufficient.
%
% * Automated time stepping based on the time-discretisation error (10 October 2026):
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
% separately. Note that the time stepping during an initial adjustment period can have a lasting effect on the solution, so excluding this
% period through StartTime does not remove its effect. In tests (Greenland, nx50, 10 years), "error-estimate" with atol=0.001 m/yr gave
% errors at the end of the run comparable to, or slightly smaller than, "iteration-based" with ATSTargetIterations=4, at about 10% lower
% computational cost. A test with a slightly different norm, corresponding to an atol of roughly 0.0006 to 0.0008 m/yr with the present
% norm, gave errors 2 to 3 times smaller than with ATSTargetIterations=4 at about the same cost. The largest remaining errors were at
% fast-flowing outlets, where local mesh refinement may also be needed.
%
% The previous parameter CtrlVar.ATSTimeDiscretisationError.Use (only used for a few days) is no longer used; if set to true,
% CtrlVar.AdaptiveTimeSteppingMethod is set to "error-estimate".
%
% *Change Log* _September 2026_
%
% * The conj grad Ua optimisation completely rewritten, and now seem to be competitive with the MATLAB lBFGS approach.
% Gradients used in optimisation are Sobolov gradients with respect to a metric defined by the explicit Hessians of the
% regularization part. Selecting inversion algorithms, e.g. gradient-based versus Hessian-based has been streamlined and
% simplified. This may possibly require some changes in input files. Box constraints in inversion can be implemented using
% two-sided variable transform. UaHessian inverse-solver uses trust-region approach with a 2D exact sub-space minimisation.
%
%
% *  UK spelling used for `Optimisation' 
%
% *Change Log* _August 2026_
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
% (No longer the case since 9 October 2026, K is again a rate. See the October 2026 Change Log.)
%
% *Change Log* _July 2026_
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
% *Change Log* _June 2026_
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
% *Change Log* _March 2026_
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
% *Change Log* _January 2026_
%
% When calculating $\dot{h}$ explicitly,  homogenized  thickness ($h$) boundary conditions are applied to the $\dot{h}$ solve. So if
% thickness is set to some prescribed value at some nodes, $\dot{h}$ at those nodes will be forced to be equal to zero. And if
% thickness links/ties are defined between nodes, those same ties will be applied to $\dot{h}$. This has no effect on $\dot{h}$
% values obtained in transient runs. But this will affect calculated/modeled $\dot{h}$ values used in an inversion when compared
% against measured $\dot{h}$ values, in combination with thickness boundary conditions.
%
%
%
% *Change Log* _December 2025_
%
% * When using $\dot{h}$ measurements in with boundary conditions applied to h, the adjoint now uses the $\dot{h}$ calculated
% explicitly using those BCs. Previously, when evaluating the sensitivities of the $\dot{h}$-cost function term at integration
% points as part of the adjoint approach, this was not the case. 
%
% 
%
% *Change Log* _November 2025_
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
% *Change Log* _July 2025_
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
% *Change Log* _May 2025_
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
% *Change Log* _March 2025_
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
% *Change Log* _September 2024_
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
% *Change Log* _July 2024_
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
% *Change Log* _March 2024_
%
% * Positive thickness constraints --- added dynamically when using the active-set option --- are now not added to nodes already
% contained in a user-defined constraint. Hence, any user-defined thickness constraints are respected, even if this means
% that some thickness go below the minimum specified thickness. Previously, thickness constraints were added to all nodes
% with thickness less than CtrlVar.ThicMin.
% 
%
%
% *Change Log* _February 2024_
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
% *Change Log* _October 2023_
%
% * A (rare) case where the Cauchy direction is not a direction of descent was incorrectly updated. This has now been addressed. 
%
%
% *Change Log* _September 2023_
%
% * Thanks to Sebastian Rosier for providing an example where the backtracking algorithm failed. This has now been corrected.
%
%
% *Change Log* _July 2023_
%
% * uv and uvh solver now uses dog-leg search if Newton back-tracking results in small steps.
%
% * Thanks to Sainan Sun for spotting that in an adaptive mesh step, call to calving was ahead of call to geometry and
% densities. This has now been corrected.
%
% *Change Log* _July 2022_
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
% *Change Log* _Sept 2021_
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
% *Change Log*
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
% *Change Log*
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
