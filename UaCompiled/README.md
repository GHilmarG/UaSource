# UaCompiled: compiled uvh assembly

Created 25 October 2026.

## Use

Set `CtrlVar.uvhMatrixAssembly.CompiledKernel=true` (default `false`). `uvhMatrixAssembly.m` then calls `uvhCompiledKernelCheck.m`, which decides for each assembly whether the compiled assembly can be used: element type with at most 10 nodes and 37 integration points, the cached assembly pattern `MUA.uvhAssemblyPattern`, the MEX files, and options supported by the kernel (see below). If so, `uvhMatrixAssemblyCompiled.m` assembles the residual and the Jacobian; if not, the vectorised assembly is used and the reason is printed once per MATLAB session. Within spmd workers (`CtrlVar.Parallel.uvhAssembly.spmd.isOn=true`) the compiled assembly is not used (`uvhMatrixAssemblySSTREAM_SPMD.m` switches it off for the workers). The results are the same as those of the vectorised assembly to round-off.

For larger meshes the fastest combination found was the compiled uvh assembly together with spmd for the uv assembly:

    CtrlVar.uvhMatrixAssembly.CompiledKernel=true;
    CtrlVar.Parallel.uvhAssembly.spmd.isOn=false;
    CtrlVar.Parallel.uvAssembly.spmd.isOn=true;

## How it works

* `uvhElementKernel.m`: loops over elements (`parfor`, OpenMP), with the loop over integration points inside, and gathers the nodal values of each element from the nodal vectors and the connectivity. Contains the terms of `uvhAssemblyIntPointImplicitSUPG_v2.m` and `uvhNodalLoopSSTREAM_v2.m`. Returns the element residuals (`Tx ... Fh`, Nele x nod) and the nonzero values of the Jacobian in CSC order: the element-matrix values are kept inside the MEX file and summed into the nonzeros in parallel over the nonzeros (using the inverse map from `uvhCSCPattern.m`). Basal drag, effective viscosity and floating mask at the integration points are only returned if requested.
* `uvhMatrixAssemblyCompiled.m`: prepares the inputs, calls the kernel, sums the residual with the precomputed indices of the assembly pattern (`P.iR`), and creates the sparse matrix directly from its CSC arrays (`SparseFromCSC`). This replaces `accumarray` and `sparse`. Entries that sum to exactly zero are kept as stored zeros.
* Data movement: `MUA.Deriv` and `MUA.DetJ` are passed to the kernel as separate arrays (fields of a structure input are copied into a MEX file, arrays are not), outputs that are completely overwritten are not zero-filled (`coder.nullcopy`), and the element-matrix values never leave the kernel MEX file.
* `uvhCSCPattern.m`: per-mesh data derived from `MUA.uvhAssemblyPattern` (inverse map for summing into the nonzeros, CSC arrays, positions of the diagonal), kept in a persistent cache keyed by the topology key of the mesh.

## Files

| File | Purpose |
|---|---|
| `uvhCompiledKernelCheck.m` | Called by `uvhMatrixAssembly.m`: decides if the compiled assembly can be used, and selects the MEX files |
| `uvhMatrixAssemblyCompiled.m` | The compiled assembly path (see above) |
| `uvhElementKernel.m` | The kernel |
| `uvhKernelFull.m`, `uvhKernelRonly.m` | Entry points for code generation: residual and Jacobian / residual only (`Ronly` and the derivative flags are compile-time constants) |
| `MakeKernelParamsFlex.m` | Run-time parameter structure from `CtrlVar` (fixed set of fields, integer codes for the sliding law and the SUPG tau option); `ok=false` if the options are not supported |
| `uvhCSCPattern.m` | Per-mesh data for the summation into the nonzeros (in the kernel) and for `SparseFromCSC` |
| `SparseFromCSC.c` | Creates a sparse matrix from its CSC arrays (hand-written C MEX: `SparseFromCSC`) |
| `BuildUaCompiled.m` | Builds the MEX files |
| `UaCompiledBuildInfo.m`, `UaCompiledBuildInfo_<mexext>.mat` | Build information: MATLAB release and checksums of the source files of each MEX file, recorded by `BuildUaCompiled.m`, checked by `uvhCompiledKernelCheck.m` once per session; outdated MEX files are not used |
| `Tests/TestUaCompiled.m` | Regression test (matlab.unittest): compiled against vectorised assembly for all options, element types and quadrature rules, fallback, one Newton solve, build information. `runtests("TestUaCompiled")`, about 25 s |
| `uvhKernel{Full,Ronly}_mex` | Kernel, any element type with nod<=10 and nip<=37 |
| `uvhKernel{Full,Ronly}_n3i6_mex` | Kernel, linear elements with 6 integration points (fixed sizes, faster) |

## Building

If the MEX files are outdated (a source file changed since they were built, or another MATLAB release), or fail when called (eg "Invalid MEX-file"), the vectorised assembly is used, with a message saying why; rebuild them with `BuildUaCompiled.m`.

Requires MATLAB Coder and a C compiler, preferably with OpenMP (tested with MinGW64 14.2, R2026b, Windows). Build in a fresh MATLAB session, because a MEX file with OpenMP that has been called cannot be replaced in the same session:

    matlab -batch "addpath(genpath('C:\cygwin64\home\Hilmar\Ua\UaSource')); BuildUaCompiled"

## Supported options

Run-time (no rebuild needed): all numeric parameters including `dt`; melange terms; `OnlyCalcBasalDragAndEffectiveViscosity`; all sliding laws of `BasalDrag.m`, with all their names and aliases (`SlidingLawCode.m`); level-set mass-balance feedback; thickness penalty (SoftPlus, polynomial, exponential); thickness barrier; `ResetThicknessToMinThickness`, `Enforce_bAboveB`; SUPG tau options tau1, tau2, taus, taut, tau0; linear, quadratic and cubic elements; quadrature rules with at most 37 points.

Not compiled (the vectorised assembly is used): quadrature rules with more than 37 points, meshes without a cached assembly pattern, the per-integration-point debugging plots (`doplots` with `PlotSUPGparameter`, `InfoLevelThickMin>=10`).

## Functions from UaMain used by the kernel

The kernel calls `BasalDrag.m` (with `SlidingLawBeta.m`), `SUPGtau.m`, `ThicknessPenaltyMassBalanceFeedback.m`, `EffectiveViscositySSTREAM.m` and other pointwise functions of UaMain directly; MATLAB Coder compiles them from there, so a change to one of them reaches the compiled kernel when the MEX files are rebuilt. They are written so that Coder can use them: options as integer codes where they are compared at every integration point (`SlidingLawCode.m`, `SUPGtauCode.m`, set once per assembly in `MakeKernelParamsFlex.m`), plotting, warnings and diagnostic output excluded from compiled code with `coder.target('MATLAB')`, no string arrays, all outputs defined on all paths. (Until 26 Oct 2026 the kernel used modified copies of `BasalDrag.m`, `SUPGtau.m` and `ThicknessPenaltyMassBalanceFeedback.m`.)

## Tests (25 October 2026)

Run the regression test after any change to the kernel, to the functions it compiles, or to the vectorised uvh assembly, and after rebuilding the MEX files:

    results=runtests("TestUaCompiled"); table(results)

It builds a small synthetic state itself (no files needed), checks for each case that the compiled path is actually used, and compares residual and Jacobian with the vectorised assembly (relative tolerance 1e-11). 40 tests, about 25 s.

With `CompiledKernel=false` the results are bit-identical to the vectorised assembly. With `CompiledKernel=true` residual and Jacobian agree with the vectorised assembly to rounding errors for all options listed above, for linear, quadratic and cubic elements, and for quadrature rules with 1 to 33 points.

* Thule example, full run (0.1 yr, 51 time steps, mesh adaptation, level set, active set): identical time steps and Newton iterations, solution identical to rounding errors (max |dh| 9e-12 m); uvh assembly 91 s -> 17 s, run time 209 s -> 146 s (profiled).
* Greenland example, one assembly, serial vectorised / spmd vectorised (4 threads) / compiled: nx200 (197,204 elements) 3903 / 1681 / 244 ms for residual and Jacobian, 1433 / 477 / 93 ms for the residual only (ratios 16.0 : 6.9 : 1 and 15.4 : 5.1 : 1).
* Greenland example, nx200, full run (0.1 yr), compiled uvh assembly with spmd for the uv assembly: 133 s, against 189 s with spmd (4 threads) for both and 326 s serial vectorised; same time steps, Newton iterations and solution to rounding errors.

Test files are in `UaExamples/Thule/PerfAnalysis`.
