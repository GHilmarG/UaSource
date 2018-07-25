function Ua(UserVar,varargin)


%% Úa
% A finite-element ice-flow model.
% 
% Flow approximations: Shallow Ice Sheet Flow Approximation (SSHEET or SIA)
%                      Shallow Ice Stream Flow Approximation (SSTREAM or SSA)
%                      Hybrid combinations of SIA and SSA (already implemented,
%                      but still being developed further...)
%
% The code is developed by Hilmar Gudmundsson (hilmar.gudmundsson@northumbria.ac.uk)
%
%% Running Úa
%
% 1) Add the folder with the Úa m-files, and its subfolders, to your matlab
% path. This can be done using the 'Home/Set' Path menu item, or from the
% command prompt doing something like:
%
%   addpath(genpath('MyUaSourceFileFolder'))
%
% 2) Define the Matlab environmental variable 'UaHomeDirectory'.
%    This can for example be done as follows:
%
%   setenv('UaHomeDirectory','MyUaSourceFileFolder')
%
% 3) If using the mesh generator `gmsh' (almost always the case) then also
%    define the Matlab environmental variable 'GmshHomeDirectory'. The gmsh
%    program for windows is in a subfolder of Ua/Source So if you are running
%    windows
%
%   setenv('GmshHomeDirectory','MyDrive/Ua/Source/gmsh-2.12.0-Windows')
%
% will do. Alternatively you might want to install your own copy of gmsh. If
% running on a Unix system, then most likely gmsh can be called without the need
% to set the Matlab environmental variable 'GmshHomeDirectory'.
%
% Now you can run Úa from within Matlab by writing:
%
%  Ua [Ret]
%
% Summary: Before running Ua do for example:
%
%   setenv('UaHomeDirectory','C:\cygwin64\home\Hilmar\ghg\Ua\Source')
%   setenv('GmshHomeDirectory','C:\cygwin64\home\Hilmar\ghg\Ua\Source\gmsh-2.12.0-Windows')
%   UaHomeDirectory=getenv('UaHomeDirectory'); addpath(genpath(UaHomeDirectory))
%
% Then run Ua from the Source directory to see if everything is OK.
%% Getting help
% You can get help on the use of Úa in the same way as you would get help on
% various in-build matlab commands by writing `help Ua'  in the matlab command
% line,  or 
%
%   doc Ua 
%
% Most m-files that are part of the Ua program have some
% inbuilt help text, for example try 
%
%   doc Ua2D_DefaultParameters
%
%
%% Defining model run
%
% Whenever setting up your own model, create your own working directory for your
% model runs. Most of the parameters of a given model run are defined by the
% user through ` *user m-files'* that are called by Úa during the run.
%
%% User m-files :
% 
% * Ua2D_InitialUserInput
% * DefineAGlenDistribution.m
% * DefineSlipperyDistribution.m
% * DefineBoundaryConditions.m 
% * DefineGeometry.m
% * DefineDensities.m
% * DefineMassBalance.m
% * DefineDesiredEleSize.m
% * DefineStartVelValues.m
% * DefineInputsForInverseRun.m
% * DefineInverseModellingVariables.m
% * UaOutputs.m
%
% To get further information on how to use individual user m-files use help. For
% example: help DefineGeometry 
% Make sure to do this from the Úa home directory, or at least not from another
% directory that has a m-file with the same name.
%
%
% When defining a new model-run, just copy these user m-files from the Úa home
% directory into your own model-run directory and modify as needed.
%
% Not all of these user m-files are always needed. For example
% `DefineInverseModellingVariables.m' is only needed for inverse runs. Also
% 'DefineStartVelValues.m' is not often required as setting start velocities to
% zero (the default opton) is usually a good approach. And DefineDesiredEleSizes
% is only needed if one finds that the standard remeshing options within Úa are
% too limited.
%
% UaOutputs is only needed for producing output files or for some plotting, etc.
%
% If any of the above listed m-Files are not found in the run directory, the
% corresponding m-Files in the Úa home directory are used instead.
%
%
%% Name of variables
% 
% Throughout the following variables stand for:
%
%  s          : upper glacier surface elevation 
%  b          : lower glacier surface elevation 
%  S          : ocean surface elevation 
%  B          : bedrock elevation 
%  (ub,vb,wb) : sliding velocity components 
%  (ud,vd,wd) : deformational  velocity components 
%  rho        : ice density (defined at nodes and does not have to be spatially uniform) 
%  rhow       : ocean density (a scalar) 
%  AGlen      : rate factor of Glen's flow law  (either a nodal or an element variable)
%  n          : stress exponent of the Glen's flow law  
%  C          : basal slipperiness ((either a nodal or an element variable) 
%  m          : stress exponent of the basal sliding law (a scalar)
%  as, ab     : mass balance distribution over the upper (as) and lower (ab) glacier surfaces. 
%               The mass balance is given in units distance/time, and should be
%               in same units as the velocity.
% alpha       : Slope of the coordinate system with respect to gravity.
% g           : The gravitational acceleration.
% GF          : a floating/grounded mask. This is a structure with the two
%                fields: node and ele. This are 1 if a node/element is grounded, 0 if
%                node/element is afloat.
%    
%
%
%
%% The variable CtrlVar
% Úa uses a the variable `CtrlVar' to define various run parameters.
%
% This variable is defined by the user in `Ua2D_InitialUserInput.m'.
% The variable has a large number of fields. List of all fields with
% descriptions can be found in  `Ua2D_DefaultParameters.m' 
%
%   help Ua2D_Defaultparameters
%   doc Ua2D_DefaultParameters
%
% CtlrVar is only defined at the start of the run in Ua2D_InitialUserInput.
% It can not be modified in any of the other user m-input files.
%
% In a restart run, CtrlVar is again always defined at the beginning of the run
% in Ua2D_InitialUserInput. It is not read from the restart file.
% 
%% The variable UserVar
%
% The variable UserVar is given as an input and output to all user m-files.
%
% It is never used by Úa. 
%
% The user can use this variable to exchange information between his own user
% m-files, or for whatever other purpose required. 
%
% UserVar can be modified at the start of the run in Ua2D_InitialUserInput.
% as well as in all user m-input files.
%
% In a restart run, UserVar is again always defined at the beginning of the run
% in Ua2D_InitialUserInput. It is not read from the restart file.
%
%% Getting information about the FE mesh from within user m-files:
%
% In all user m-files the variable MUA is given as an input. MUA contains all
% information about the FE mesh.
%
% MUA is a structured variable with the following fields:
%
%      coordinates:  Nnodes x 2 array with the x and y coordinates of all nodal
%                    points
%     connectivity:  mesh connectivity
%           Nnodes:  number of nodes in mesh
%             Nele:  number of elements in mesh
%              nod:  number of nodes per element nip:  number of integration
%                    points
%           points:  local element coordinates of integration points
%          weights:  weights of integration points
%         Boundary:  a structure containing info about mesh boundary
%                    This structure is calculated as:
%                    MUA.Boundary=FindBoundary(connectivity,coordinates); and
%                    info about the fields can be found in `FindBoundary.m'
%            Deriv:  element derivatives
%             DetJ:  element determinants
%
% The values of MUA should never be changed directly by the user.
%
%
%% Meshing
% There are various ways of meshing the computational domain.
%
% In almost all cases the simplest option tends to be to define the outlines of
% the computational domain in Ua2D_InitialUserInput. In that case Úa will call
% an external mesh generator. The external mesh generator used by Úa is "gmsh"
% which is a well known and well supported open source mesh generator
% (http://geuz.org/gmsh/). The outlines of the mesh are defined by the variable
% 'MeshBoundaryCoordinates' set in Ua2D_InitialUserInput.m. This approach is
% quite flexible and allows for complicated computational domains containing
% holes and/or separated domains. 
%
% For examples of how to generate different
% type of meshes run *ExamplesOfMeshGeneration*
%
% The ExamplesOfMeshGeneration.m contains information and examples on how
% to define inputs for various types of meshes. 
%
% There are also various ways of refining the mesh. Both global and local
% (explicit) adaptive meshing is supported. See further explanations in
% 'Ua2D_DefaultParamters.m'
%%

%
%


if nargin==0
    UserVar=[];
end

Ua2D(UserVar,varargin{:})


end


