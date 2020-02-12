



%%
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
