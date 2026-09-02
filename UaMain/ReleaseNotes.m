

%%
%
% *Release Notes* _August 2026_
%
%
% * For the -uvh- solve, two normalization of the residuals are now available, selected by
%
%    CtrlVar.uvhResidualNormalisation = "pooled"    (default, historical)
%                                     = "blockwise"
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
% * MATLAB seems to have been busy working on their optimization functions, and the performance of using
%
%    CtrlVar.Inverse.MinimisationMethod="MatlabOptimization-GradientBased";  
%
% now appears improved. This option actually stop working in Matlab 2021b and 2022a, and this may have been due to
% a bug in the optimisation toolbox. From at least 2024a onward this now works again, and based on some numerical tests,
% appears much improved. This is currently not the default option, but users might consider setting
%
%    CtrlVar.Inverse.MinimisationMethod="MatlabOptimization-GradientBased";  
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
