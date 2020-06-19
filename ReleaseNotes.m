



%%
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
% You can also now specify a minimum selected time step in the
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
