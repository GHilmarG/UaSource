
function Ua




%% Úa
% A finite-element ice-flow model.
% 
% Flow approximations: Shallow Ice Sheet Flow Approximation (SSHEET or SIA)
%                      Shallow Ice Stream Flow Approximatoin (SSTREAM or SSA)
%                      Hybrid combinations of SIA and SSA (already implemented, but still being developed further...)
%
%% Running Úa:
%
% 1) Add the folder with the Úa m-files, and its subfolders, to your matlab path.
%    This can be done using the 'Home/Set' Path menu item, or from the
%    command prompt doing something like:
%    addpath('MyUaSourceFileFolder')
%
% 2) Define the Matlab environmental variable 'UaHomeDirectory'.
%    This can for example be done as follows:
%                 setenv('UaHomeDirectory','MyUaSourceFileFolder')
%
% 3) If using the mesh generator `gmsh' (almost always the case) then also define 
%    the Matlab environmental variable 'GmeshHomeDirectory'. The gmsh program
%    can be found below the Úa main folder, so something like
%                 setenv('GmeshHomeDirectory','MyDrive/Ua2D/gmsh-2.8.4-Windows')  
%   will do. Alternativily you might want to install your own copy of gmsh.
%   If running on a Unix system, then most likely gmsh can be called without the need
%   to set the Matlab environmental variable 'GmeshHomeDirectory'. 
%
% Now you can run Úa from within Matlab by writing Ua2D [Ret]
%
% Example:
%
%   setenv('UaHomeDirectory','C:\cygwin64\home\Hilmar\ghg\Ua\Source\Ua2D_Development') 
%   setenv('GmeshHomeDirectory','C:\cygwin64\home\Hilmar\ghg\Ua\gmsh-2.11.0-Windows')
%   UaHomeDirectory=getenv('UaHomeDirectory');
%   addpath(genpath(UaHomeDirectory))
%
%
%% Getting help
% You can get help on the use of Úa in the same way as you would get help on various in-build matlab commands by writing 
% `help Ua'  in the matlab command line,  or `doc Ua'.
% Most m-files that are part of the Ua program have some inbuild help text, for example try `doc Ua2D_DefaultParameters'.
%
% To get a HTML formatted documentation try:
%
%    publish('Ua.m','evalCode',false) ; web('html\Ua.html')
%
%% Defining model run
%
% Whenever setting up your own model, create your own working directory for your model runs. 
% Most of the parameters of a given model run are defined by the user through ` *user m-files'* 
% that are called by Úa during the run. 
%
%% User m-files :
% 
% * Ua2D_InitialUserInput
% * DefineAGlenDistribution.m
% * DefineSlipperyDistribution.m
% * DefineBoundaryConditions.m  (also possible to use the more limited but easier to use `DefineBCs.m' instead)          
% * DefineGeometry.m                
% * DefineDensities.m          
% * DefineMassBalance.m
% * DefineDesiredEleSize.m
% * DefineStartVelValues.m
% * DefineInverseModellingVariables.m
% * UaOutputs.m
%
% 
% 
% When defining a new model-run, just copy these files from the Úa home directory into your own model-run directory 
% and modify as needed. 
%
% Not all of these user m-files are always needed. For example `DefineInverseModellingVariables.m' is only needed for inverse runs.
% Also  'DefineStartVelValues.m' is not often required as setting start velocities to zero (the default opton) is usually a good approach.
% And DefineDesiredEleSizes is only needed if one finds that the standard remeshing options within Úa are too limited.
%
% UaOutputs is only needed for producing output files or for some plotting, etc.  
%
% If any of the above listed m-Files are not found in the run directory, the corresponding m-Files in the Úa home directory are used instead.
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
%  (ud,vd,wd) : deformational velocity components
%  rho        : ice density (defined at nodes and does not have to be spatially uniform)
%  rhow       : ocean density (a scalar)
%  AGlen      : the rate factor in Glen's flow law  (either a nodal or an element variable)
%  n          : stress exponent in Glen's flow law  (either a nodal or an element variable)
%  C          : basal slipperiness (a scalar)
%  m          : stress exponent in basal sliding law (a scalar)
%
%% Calls to user m-files:
% The calls to these functions are: 
%
%    [rho,rhow,g]=DefineDensities(Experiment,CtrlVar,MUA,time,s,b,h,S,B);
% 
%    [C,m]=DefineSlipperyDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
% 
%    [AGlen,n]=DefineAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
% 
%    [as,ab]=DefineMassBalance(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF); 
%
% If mass-balance geometry feedback is included, define mass balance as:
% 
%    [as,ab,dasdh,dabdh]=DefineMassBalance(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF); 
% 
%    BCs=DefineBoundaryConditions(Experiment,CtrlVar,MUA,BCs,time,s,b,h,S,B,ub,vb,ud,vd,GF)
%
% BCs is an instant of the class `BoundaryConditions'.
% To see the fields of BCs have a look at BoundaryConditions.m 
% The fancy way of doing this is to do: publish('BoundaryConditions.m') ; web('html\BoundaryConditions.html')
% from the Úa home directory.
% Alternativily, instead of `DefineBoundaryConditioins' one can use:
%
%   [ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,hfixednode,hfixedvalue,htiedA,htiedB]=...
%     DefineBCs(Experiment,CtrlVar,MUA,time,s,b,h,S,B,ub,vb,ud,vd,GF);
% 
%   [s,b,S,B,alpha]=DefineGeometry(Experiment,CtrlVar,MUA,time,FieldsToBeDefined);
% 
%   [ub,vb,ud,vd]=DefineStartVelValues(Experiment,CtrlVar,MUA,ub,vb,ud,vd,time,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m);
%  
%

%% Getting information about mesh from within User m-files:
%
% In all user m-files the variable MUA is given as an input.
%  MUA is a structured variable with the following fields
%      coordinates:  Nnodes x 2 array with the x and y coordinates of all nodal points
%     connectivity:  mesh connectivity 
%           Nnodes:  number of nodes in mesh
%             Nele:  number of elements in mesh
%              nod:  number of nodes per element
%              nip:  number of integration points
%           points:  local element coordinates of integration points
%          weights:  weights of integration points
%         Boundary: a structure containing info about mesh boundary
%                   This structure is calculated as:
%                   MUA.Boundary=FindBoundary(connectivity,coordinates);
%                   and info about the fields can be found in `FindBoundary.m'
%            Deriv:  element derivatives
%             DetJ:  element determinants
%
%               The values of MUA should never be changed directly by the user.
%%
% see also Ua2D_DefaultParameters
%%





if nargin==0
    UserRunParameters=[];
end

Ua2D(UserRunParameters)


end


