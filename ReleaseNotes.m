

%%
%
% *Release Notes* _March 2025_
%
% The uv assembly was changed slightly to make the assembly with respect to h more consistent with the uvh assembly. This can
% lead to uv solve to be different if the ice density is varying greatly specially. The basic difference is that the
% flotation condition, which involves the product rho and h, is now evaluated at integration points directly. Previously in
% the uv solve it was evaluated by first forming the product, rho h, at nodes and then interpolated to the integration
% points.
%
% Calls to gcp have been replaced with gcp('nocreate'). The implication is that a parallel pool should be defined and started
% ahead of a call to Ua. However, not that this might also depend on user settings for the parallel pool. For example if the
% user setting imply automated start of a parallel pool whenever parfor or smpd is encountered. If no parallel pool is found,
% all parallel options are turned off, and the run proceeds in non-parallel mode.  
%
% Unless the mass balance/thickness feedback is activated, MassBalance evaluation now lags behind by one time step. This was
% done to reduce the calls to DefineMassBalance, and to make sure that the mass balance of the Ua fields, F0, was same as the
% F from previous time step.  Note that when the mass balance/thickness feedback is activated, the mass balance function is
% called within the assembly loop and then the mass balance will not lag behind. 
%
% The start and end times can now be specified using (new option):
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
% When not inverting for both A and C, the variable InvValues was not updated in last call, leading to a possible mismatch
% between F.A and InvValues.A. Thanks to Camilla Schelpe for identifying this.
%
% The implementation of the ice-thickness barrier function has been simplified, and is not similar to the barrier function
% used in the level-set solver.
%
% The MassContinuity solver, (only used when solving for h alone, and not in uv or uvh solves), not uses the active-set
% method to enforce positive thickness.
%
% The model was tested with MATLAB R2024b 
%
% *Release Notes* _September 2024_
%
% For comparison purposes the semi-implicit solver has been updated. 
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
% The default value of 
%
%   CtrlVar.LevelSetMethodThicknessConstraints;
%
% is now set to true (previously set to false). Therefore if one uses the active-set method and the level-set method, min
% thickness constraints will be automatically applied to all nodes downstream of calving fronts, unless of course this
% parameter is set to false by the user. 
% 
% *Release Notes* _July 2024_
%
% The utility
%   
%   [tbx,tby,tb,eta] = CalcBasalTraction(CtrlVar,UserVar,MUA,F,options)
%
% has been updated to allow for calculation at integration points (before only calculated nodal values).
%
% An inconsistent input parameter test when using different sliding laws for basal drag and ocean drag calculations
% corrected. Thanks to Sainan Sun for point this out. 
%
% A situation where a composite variable was used inside of spmd, resulting in a warning but no errors, has been corrected.
%
% *Release Notes* _March 2024_
%
% The utility
%   
%   [tbx,tby,tb,eta] = CalcBasalTraction(CtrlVar,UserVar,MUA,F,options)
%
% has been updated to allow for calculation at integration points (before only calculated nodal values).
%
% An inconsistent input parameter test when using different sliding laws for basal drag and ocean drag calculations
% corrected. Thanks to Sainan Sun for point this out. 
%
% A situation where a composite variable was used inside of spmd, resulting in a warning but no errors, has been corrected.
%
% *Release Notes* _March 2024_
%
% Positive thickness constraints --- added dynamically when using the active-set option --- are not added to nodes already
% contained in a user-defined constraint. Hence, any user-defined thickness constraints are respected, even if this means
% that some thickness go below the minimum specified thickness. Previously, thickness constraints were added to all nodes
% with thickness less than CtrlVar.ThicMin.
%
%
%
% *Release Notes* _February 2024_
%
% Parallel spmd assembly option now improved and shows good scalability, although speedup always somewhat problem dependent.  
% On local workstations with 12 workers, speedup ranging from 6 to 10 seems easily obtainable, and on machines with 24 workers a
% speedup of 20 has been observed. 
%
% To switch on the smpd parallel assembly do:
%
%   CtrlVar.Parallel.uvhAssembly.spmd.isOn=true;    % uvh assembly in parallel using spmd over sub-domain (domain decomposition)  
%   CtrlVar.Parallel.uvAssembly.spmd.isOn=true;     % uv assembly in parallel using spmd over sub-domain (domain decomposition)  
%
%
% The linear system can now also be solved using distributed arrays, although not all cases yet implemented. Turn this on/off as:
%
%   CtrlVar.Parallel.Distribute=true/false;                       % linear system is solved using distributed arrays. 
%
% Again, the speedup is problem dependent and never particularly large for low density sparse matrices as typically generated
% by Ua and other FE programs. 
%
%   CtrlVar.Parallel.isTest=false/true;             % Runs both with and without parallel approach, and prints out some information on relative performance. 
%                                                   % Good for testing if switching on the parallel options speeds things up, and by how much.
% 
% *Release Notes* _October 2023_
%
% A (rare) case where the Cauchy direction is not a direction of descent was incorrectly updated. This has now been addressed. 
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
% Cornford, and  Umbi. Refer to the Ua compendium for definition of these
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
