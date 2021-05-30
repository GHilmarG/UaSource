



%%
%
% *Release Notes* _June 2021_
%
% * Default inverse algorithim has changed. Now the defaul options is a Hessian-based inversion.  The older-approach is still available
% by selecting a Gradien-based inversion. This is specificed using the CtrlVar field:
%
%       CtrlVar.Inverse.MinimisationMethod
%
% When using this Hessian-based inversion, there is no need to specify the value of 
%
%       CtrlVar.Inverse.AdjointGradientPreMultiplier 
%
% and this field is not used.
%   
% * Ocean and wind-indued drag can be included. This is defined in 
%
%   DefineSeaIceParameters.m
%
% and you need to switch this option on by setting
%
%   CtrlVar.IncludeMelangeModelPhysics=true;
%
% in DefineInitialInputs.m
%
% * The quadrature degree can now be directly specified up to degree 25.  You specify this using
%
%   CtrlVar.QuadratureRuleDegree
%
% Generally, the default option should be fine, but there might be instances were increasing the degree could help with obtaining
% second-order convergence of the Newton-Raphson system.
%
% * The shallow-ice sheet (SIA/SSTREAM) option now includes basal sliding.
%
% The SSHEET option uses the Weertman sliding law only. The transient SSHEET solution is done implicitly with respect to the thickneess
% using the NR method.  When using SSTREAM you will, in general, need to specify boundary conditions for both the deformational and the
% basal sliding velocities.
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
% If you define these old fields in your (new) DefineInitialInputs.m, Úa will spot this and complain bitterly. 
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
% * Inversion can be done using all implemented sliding laws except Coulumb. 
% 
%
% * u, v and h residuals now calulated in the L2 norm instead of the l2 norm as in the
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
