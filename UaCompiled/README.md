# UaCompiled: compiled uvh assembly

Created 25 October 2026.

## Use

Set `CtrlVar.uvhMatrixAssembly.CompiledKernel=true` (default `false`). `uvhMatrixAssembly.m` then calls `uvhCompiledKernelCheck.m`, which decides for each assembly whether the compiled assembly can be used: element type with at most 10 nodes and 37 integration points, the cached assembly pattern `MUA.uvhAssemblyPattern`, the MEX files, and options supported by the kernel (see below). If so, `uvhMatrixAssemblyCompiled.m` assembles the residual and the Jacobian; if not, the vectorised assembly is used and the reason is printed once per MATLAB session. Within spmd workers (`CtrlVar.Parallel.uvhAssembly.spmd.isOn=true`) the compiled assembly is not used (`uvhMatrixAssemblySSTREAM_SPMD.m` switches it off for the workers). The results are the same as those of the vectorised assembly to round-off.

For larger meshes the fastest combination found was the compiled uvh assembly together with spmd for the uv assembly:

    CtrlVar.uvhMatrixAssembly.CompiledKernel=true;
    CtrlVar.Parallel.uvhAssembly.spmd.isOn=false;
    CtrlVar.Parallel.uvAssembly.spmd.isOn=true;

## How it works

* `uvhElementKernel.m`: loops over elements (`parfor`, OpenMP), with the loop over integration points inside, and gathers the nodal values of each element from the nodal vectors and the connectivity. Contains the terms of `uvhAssemblyIntPointImplicitSUPG_v2.m` and `uvhNodalLoopSSTREAM_v2.m`. Returns the element residuals (`Tx ... Fh`, Nele x nod), the element-matrix values element-major (`Xmat`, 9*nod^2 x Nele), and basal drag, effective viscosity and floating mask at the integration points.
* `uvhMatrixAssemblyCompiled.m`: prepares the inputs, calls the kernel, sums the residual with the precomputed indices of the assembly pattern (`P.iR`), sums the element-matrix values directly into the nonzeros of the sparse matrix (`GatherSumCSC`) and creates the sparse matrix directly from its CSC arrays (`SparseFromCSC`). This replaces `accumarray` and `sparse`. Entries that sum to exactly zero are kept as stored zeros.
* `uvhCSCPattern.m`: per-mesh data derived from `MUA.uvhAssemblyPattern` (inverse map, CSC arrays, positions of the diagonal), kept in a persistent cache keyed by the topology key of the mesh.

## Files

| File | Purpose |
|---|---|
| `uvhCompiledKernelCheck.m` | Called by `uvhMatrixAssembly.m`: decides if the compiled assembly can be used, and selects the MEX files |
| `uvhMatrixAssemblyCompiled.m` | The compiled assembly path (see above) |
| `uvhElementKernel.m` | The kernel |
| `uvhKernelFull.m`, `uvhKernelRonly.m` | Entry points for code generation: residual and Jacobian / residual only (`Ronly` and the derivative flags are compile-time constants) |
| `MakeKernelParamsFlex.m` | Run-time parameter structure from `CtrlVar` (fixed set of fields, integer codes for the sliding law and the SUPG tau option); `ok=false` if the options are not supported |
| `uvhCSCPattern.m` | Per-mesh data for `GatherSumCSC` and `SparseFromCSC` |
| `GatherSumCSC.m` | Sums the element-matrix values into the nonzeros (compiled: `GatherSumCSC_mex`) |
| `SparseFromCSC.c` | Creates a sparse matrix from its CSC arrays (hand-written C MEX: `SparseFromCSC`) |
| `BasalDragFlex.m`, `SUPGtauFlex.m`, `ThicknessPenaltyMassBalanceFeedbackFlex.m` | Modified copies of the UaMain functions, see below |
| `BuildUaCompiled.m` | Builds the MEX files |
| `uvhKernel{Full,Ronly}_mex` | Kernel, any element type with nod<=10 and nip<=37 |
| `uvhKernel{Full,Ronly}_n3i6_mex` | Kernel, linear elements with 6 integration points (fixed sizes, faster) |

## Building

Requires MATLAB Coder and a C compiler, preferably with OpenMP (tested with MinGW64 14.2, R2026b, Windows). Build in a fresh MATLAB session, because a MEX file with OpenMP that has been called cannot be replaced in the same session:

    matlab -batch "addpath(genpath('C:\cygwin64\home\Hilmar\Ua\UaSource')); BuildUaCompiled"

## Supported options

Run-time (no rebuild needed): all numeric parameters including `dt`; melange terms; `OnlyCalcBasalDragAndEffectiveViscosity`; the sliding laws of the forward branch of `BasalDrag` (Weertman/W, Budd/W-N0, Coulomb/C, Tsai/minCW-N0, Cornford/rpCW-N0, Umbi/rCW-N0, Joughin/rCW-v0); level-set mass-balance feedback; thickness penalty (SoftPlus, polynomial, exponential); thickness barrier; `ResetThicknessToMinThickness`, `Enforce_bAboveB`; SUPG tau options tau1, tau2, taus, taut, tau0; linear, quadratic and cubic elements; quadrature rules with at most 37 points.

Not compiled (the vectorised assembly is used): the inverse-run branch of `BasalDrag` (`CtrlVar.Inverse.dFuvdClambda`), the sliding-law aliases "B" and "-C-" (which the forward branch of `BasalDrag` does not accept either), quadrature rules with more than 37 points, meshes without a cached assembly pattern, the per-integration-point debugging plots (`doplots` with `PlotSUPGparameter`, `InfoLevelThickMin>=10`).

## Changes in the copies of the UaMain functions

* `BasalDragFlex.m`: `ismember(string(...),[...])` replaced by an integer code (code generation does not support string arrays, and comparing strings at every integration point was slow); `t=[]; beta2i=[]; Dbeta2i=[]` replaced by zero arrays (same size in all branches); the inverse branch is excluded with `coder.target('MATLAB')`; the forward `switch` uses the integer code `CtrlVar.SlidingLawCode`.
* `SUPGtauFlex.m`: plotting excluded with `coder.target('MATLAB')`; the `switch` uses a numeric code for the tau option.
* `ThicknessPenaltyMassBalanceFeedbackFlex.m`: plotting/info output and the warning that uses a persistent variable excluded with `coder.target('MATLAB')`.

## Tests (25 October 2026)

With `CompiledKernel=false` the results are bit-identical to the vectorised assembly. With `CompiledKernel=true` residual and Jacobian agree with the vectorised assembly to rounding errors for all options listed above, for linear, quadratic and cubic elements, and for quadrature rules with 1 to 33 points.

* Thule example, full run (0.1 yr, 51 time steps, mesh adaptation, level set, active set): identical time steps and Newton iterations, solution identical to rounding errors (max |dh| 9e-12 m); uvh assembly 91 s -> 24 s, run time 209 s -> 152 s (profiled).
* Greenland example, one assembly, serial vectorised / spmd vectorised (4 threads) / compiled: nx200 (197,204 elements) 3894 / 1628 / 386 ms for residual and Jacobian, 1417 / 490 / 159 ms for the residual only.

Test files are in `UaExamples/Thule/PerfAnalysis`.
