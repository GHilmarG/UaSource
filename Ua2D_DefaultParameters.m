function CtrlVar=Ua2D_DefaultParameters


%% Running Úa:
%
% 1) Add the folder with the Úa m-files, and its subfolders, to your matlab path.
%    This can be done using the 'Home/Set' Path menu item, or from the
%    command prompt doing something like:
%    addpath('MyUaSourceFileFolder')
%
% 2) Define the Matlab environmental variable 'UaHomeDirectory'.
% This can for example be done as follows:
%                 setenv('UaHomeDirectory','MyUaSourceFileFolder')
%
% 3) If using the mesh generator `gmsh' (almost always the case) then also define 
% the Matlab environmental variable 'GmeshHomeDirectory'. The gmsh program
% can be found below the Úa main folder, so something like
%                 setenv('GmeshHomeDirectory','MyDrive/Ua2D/gmsh-2.8.4-Windows')  
% will do. Alternativily you might want to install your own copy of gmsh.
% If running on a Unix system, then most likely gmsh can be called without the need
% to set the Matlab environmental variable 'GmeshHomeDirectory'. 
%
% Now you can run Úa from within Matlab by writing Ua2D [Ret]
%
%
%% Defining model run
% Whenever setting up your own model, create your own working directory for your model runs. 
% Most of the parameters of a given model run are defined by the user through m-files 
% called by Úa during the run. 
%
% These m-files are:
% DefineAGlenDistribution.m
% DefineSlipperyDistribution.m
% DefineBoundaryConditions.m  (also possible to use the more limited but easier to use `DefineBCs.m' instead)          
% DefineGeometry.m                
% DefineDensities.m          
% DefineMassBalance.m
%
% Optionally one can also define start values for (u,v) using:
% DefineStartVelValues.m
%
% If adaptive meshing is used, then optionally one can define desired ele sizes using 'DefineDesiredEleSize.m'
% However, Úa also has various automated in-built options of specifying ele sizes based on various criteria and
% in many cases one does not need to do any further modifications.
%
% If used in an inverse mode then 'DefineInverseModellingVariables.m' is required.
%
% If any of the above listed m-Files are not found in the run directory, the corresponding m-Files in the Úa home directory are used instead.
%
% When defining a new model-run, just copy these files from the Úa home directory into your own model-run directory 
% and modify as needed. 
% 
% The calls to these functions are: 
% [rho,rhow,g]=DefineDensities(Experiment,CtrlVar,MUA,time,s,b,h,S,B);
% 
% [C,m]=DefineSlipperyDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
% 
% [AGlen,n]=DefineAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
% 
% [as,ab]=DefineMassBalance(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF);
% or if mass-balance geometry feedback option 3 is used:
% [as,ab,dasdh,dabdh]=DefineMassBalance(Experiment,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF); 
% 
% BCs=DefineBoundaryConditions(Experiment,CtrlVar,MUA,BCs,time,s,b,h,S,B,ub,vb,ud,vd,GF)
% BCs is an instant of the class `BoundaryConditions' which in fact is just a simple structure with a number of fields. 
% To see the fields of BCs have a look at BoundaryConditions.m 
% The fancy way of doing this is to do: publish('BoundaryConditions.m') ; web('html\BoundaryConditions.html')
% from the Ua home directory.
%
% Alternativily, instead of `DefineBoundaryConditioins' one can use:
% [ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,hfixednode,hfixedvalue,htiedA,htiedB]=...
%     DefineBCs(Experiment,CtrlVar,MUA,time,s,b,h,S,B,ub,vb,ud,vd,GF);
% 
% [s,b,S,B,alpha]=DefineGeometry(Experiment,CtrlVar,MUA,time,FieldsToBeDefined);
% 
% [ub,vb,ud,vd]=DefineStartVelValues(Experiment,CtrlVar,MUA,ub,vb,ud,vd,time,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m);
%  
% [ub,vb,ud,vd]=StartVelocity(CtrlVar,MUA);
%

%% Name of variables
%  
% Throughout the following variables stand for:
%
%   s          : upper glacier surface elevation 
%   b          : lower glacier surface elevation 
%   S          : ocean surface elevation
%   B          : bedrock elevation 
%   (ub,vb,wb) : sliding velocity components
%   (ud,vd,wd) : deformational velocity components
%    rho       : ice density (defined at nodes and does not have to be spatially uniform)
%    rhow      : ocean density (a scalar)
%    AGlen     : the rate factor in Glen's flow law  (either a nodal or an element variable)
%    n         : stress exponent in Glen's flow law  (either a nodal or an element variable)
%    C         : basal slipperiness (a scalar)
%    m         : stress exponent in basal sliding law (a scalar)
%
%  
%  MUA :        The finite-element mesh structure 
%               The values of MUA should never be changed directly by the user.
%               MUA has the following fields              
%      coordinates:  Nnodes x 2 array with the x and y coordinates of all nodal points
%     connectivity:  mesh connectivity 
%           Nnodes:  number of nodes in mesh
%             Nele:  number of elements in mesh
%              nod:  number of nodes per element
%              nip:  number of integration points
%           points:  integration points
%          weights:  weights of integration points
%         Boundary: a structure containing info about mesh boundary
%                   This structure is calculated as:
%                   MUA.Boundary=FindBoundary(connectivity,coordinates);
%                   and info about the fields can be found in `FindBoundary.m'
%            Deriv:  element derivatives
%             DetJ:  elemnt determinants
%
%%

CtrlVar.fidlog=1;
CtrlVar.Experiment='UaDefaultRun';
CtrlVar.time=NaN;
%% Types of run
% 
CtrlVar.TimeDependentRun=0 ;  % any of [0|1].  
                              % If true (i.e. set to 1) then the run is a forward transient one, if not
                              % then velocities based on the current geometry are calculated. 
CtrlVar.InverseRun=0;         % if true then a surface-to-bed inversion is to be performed.
                              % (in an inverse run the value of CtrlVar.TimeDependentRun is irrelevant)

%% Ice flow approximation
CtrlVar.FlowApproximation='SSTREAM' ;  % any off ['SSTREAM'|'SSHEET'|'Hybrid']

%% Boundary conditions
CtrlVar.UpdateBoundaryConditionsAtEachTimeStep=0;  % if true, `DefineBCs' is called at the beginning of each time step and boundary conditions are updated
                                                   % otherwise boundary conditions are only updated at the beginning of the run (also at the beginning or a restart run).
CtrlVar.BCsWeights=1;  % testing parameter, do not change
% Boundary conditions are defined by the user using the m-File: DefineBCs
% if one has to define Dirichlet BCs along a complex boundary then the m-File:
%                 [Ind,AlongDist,NormDist] = DistanceToLineSegment(p, A, B,tolerance)
% can be usefull to call within DefineBC.  See comments in DistanceToLineSegment.m for explanation on how to use.
%
%%
CtrlVar.DefineOceanSurfaceAtEachTimeStep=0;   % if true,  `DefineGeometry.m' is called at each time step, returning S, and only S.
                                              % if false, `DefineGeometry.m' is only called at the beginning of a run
                                              %            and when the FE-mesh changes

%%
CtrlVar.TestUserInputs=1;  % By default user inputs will be tested at the start of the run
                           % to suppress set TestUserInputs=0
                           % if user inputs are always to be tested throughout the run, set TestUserInputs=2 
%% Element type
%
% The options are: linear, quadratic, or cubic Lagrangian triangle elements
CtrlVar.TriNodes=6 ;  % Possible values are 3, 6, 10 node (linear/quadradic/cubic)

%% Control on transient runs
% Once either the number of time steps or total time modelled reaches prescribed values
% the run stops.
CtrlVar.nTimeSteps=1;            % maximum number of time steps
CtrlVar.TotalTime=1e10;          % maximum model time
CtrlVar.dt=1;                    % time step (usually overwritten by user by defining dt in the Ua2D_InitialUserInputFile
CtrlVar.dtmin=1e-12;             % for numerical reasons the time step should always be larger than some very small value

CtrlVar.InitialDiagnosticStep=1; % start implicit transient (prognostic) run with an initial diagnostic step (a good idea, do this always)

%% Restart option
CtrlVar.Restart=0;                       % either 0/false or 1/true.  Set to 1 for a restart run
CtrlVar.WriteRestartFile=1;              % if true, a restart file is written
CtrlVar.WriteRestartFileInterval=100;    % restart file written at this time-step interval  (these are time steps, not model time)
CtrlVar.ResetTime=0 ;                    % set to 1 to reset (model) time at start of restart run
CtrlVar.RestartTime=NaN;                 % if ResetTime is true, then this is the model time at the start of the restart run
CtrlVar.ResetTimeStep=0;                 % 1 if time step should be reset to dt given in the Ua2D_InitialUserInputFile
CtrlVar.NameOfRestartFiletoRead='Ua2D_Restartfile.mat';
CtrlVar.NameOfRestartFiletoWrite='Ua2D_Restartfile.mat';

%%
CtrlVar.SaveAdaptMeshFileName=[];          % file name for saving adapt mesh. If left empty, no file is written

%% Plotting
%
% Most plotting is typically done by the user using his own version of the `UaOutputs.m',
% or in a separated post-processing step
% However, some basic plots can be generated directly from within Ua.
%

CtrlVar.doplots=1;          % if true then plotting during runs by Ua are allowed, set to 0 to suppress all plots
CtrlVar.PlotWaitBar=1;      % a waitbar is plotted
CtrlVar.doAdaptMeshPlots=1; % if true and if CtrlVar.doplots true also, then do some extra plotting related to adapt meshing
CtrlVar.TransientPlotDt=NaN;% model time interval between calls to `FE2dTransientPlots.m'
                            % special values:  
                            %              0 : always call `FE2dTransientPlots.m', i.e. at each and every time step
                            %            NaN : never call `FE2dTransientPlots.m'
                            %  Note: when using automated time stepping, the time step is adjusted to fit this criteria
CtrlVar.OnlyDoFirstTransientPlotAndThenStop=0; % stops run after first transient plot has been made, useful for initial testing purposes, this will
                                               % allow inspection of s,b,B,S,mesh, etc, before any model calculations are started

CtrlVar.PlotStrains=0;
CtrlVar.PlotOceanLakeNodes=0;        % Shows which nodes are considered a part of the `ocean' and which are within `lakes' that have no connection the ocean
CtrlVar.PlotMeltNodes=0;
CtrlVar.FE2dPlots=0;                 % if true/1 then FE2dPlots.m is called at the end of run. User would generally make his own version of this m-file
CtrlVar.PlotBackgroundImage=0;

CtrlVar.PlotXYscale=1;     % used to scale x and y axis of some of the figures, only used for plotting purposes
                           % (if spatial units are in meters, setting this to 1000 produces xy axis with the units km)
CtrlVar.PlotsXaxisLabel='x' ; CtrlVar.PlotsYaxisLabel='y' ; %

CtrlVar.PlotMesh=0;        % If true then FE mesh is shown every time a new mesh is generated
CtrlVar.FEmeshPlotTitle=[]; % Title for FE mesh plot, if left empty then something sensible is used instead
CtrlVar.PlotFEmeshAndSaveMesh=0 ; % when plotting mesh also save mesh to a file
CtrlVar.PlotBCs=0;         % If true then boundary conditions are shown at the beginning of the run
CtrlVar.BoundaryConditionsFixedNodeArrowScale=1;
CtrlVar.PlotNodes=0;       % If true then nodes are plotted when FE mesh is shown
CtrlVar.PlotLabels=0 ;     % If true elements and nodes are labelled with their respective numbers
CtrlVar.PlotNodesSymbol='o';
CtrlVar.PlotNodesSymbolSize=3;
CtrlVar.MeshColor='k'; CtrlVar.NodeColor='k';
CtrlVar.GLresolutionWhenPlotting=2000;      % when plotting GL the GF mask is (sometimes) mapped on a regular grid
CtrlVar.MinSpeedWhenPlottingVelArrows=0;    % when plotting vel arrows with smaller speed are scaled so that their speed its
                                            % equal to this value  (setting this to a large value makes all arrows
                                            % equally long)


%% Numerical variables related to transient runs
% In general there should be no need to ever change these values except for testing purposes
%
% Transient runs can be done either (fully) implicitly, or semi-implicitly
% In a (fully) implicit approach, the time-integration is done implicitly with respect to both velocities and thickness.
% In a semi-implict approach, the time-integration is done implicity with respect to thickness, and explicitly with respect to velocities.
%
% There are currently two fully-implicit time-stepping methods implemented: The 'theta' and the 'supg' methods.
%
% The 'theta' method uses a weighted sum of the values at the beginning and the end of a time step.
% The weighting is controled by CtrlVar.theta and depending on the value of theta different types of
% approximations are obtained: 0,1/2,1 gives forward Euler, Lax-Wendroff and backwards Euler, respectivily.
% The 'supg' method is a Streamline-Upwind Petrov-Galerkin method. The supg-method uses the same
% weighting as the 'theta' method, but the test function for the mass-conservation equation is different.
%
% The default time-stepping method is: Fuly implicit Streamline-Upwind Petrov-Galerkin with theta=0.5 (Lax Wendroff).
%
%
CtrlVar.Implicituvh=1;           % 0: prognostic run is semi-implicit (implicit with respect to h only)
                                 % 1: prognostic run is fully-implicit (implicit with respect to uvh)

CtrlVar.uvhTimeSteppingMethod='supg'; % 'theta'|'supg'

CtrlVar.SUPG.beta0=0.5 ; CtrlVar.SUPG.beta1=0 ; % parameters related to the SUPG method.
CtrlVar.theta=0.5;    % theta=0 is forward Euler, theta=1 is backward Euler, theta=1/2 is Lax-Wendroff and is most accurate

% Note: An additional time-stepping mehtod is the Third-Order Taylor-Galerkin (TG3) method.
% It has not been fully tested but seems to work very well for fully implicit transient calculation.
% This option that can be obtained by setting:
% CtrlVar.TG3=1 ;  CtrlVar.Test1=1;  CtrlVar.Test0=0;   CtrlVar.theta=0.5;  
% and using the fully-implicit time-stepping option (CtrlVar.Implicituvh=1)); 
CtrlVar.TG3=0 ; % if true, the prognostic steps uses a third-order Taylor-Galerkin method
                % currently only implemented for periodic boundary conditions                         
                % Note, only theta=0.5 is strictly consistent with TG3=1, so
                % for backward Euler set theta=1 and TG3=0                 
CtrlVar.IncludeTG3uvhBoundaryTerm=0;                     % keep zero (only used for testing)
CtrlVar.IncludeDirichletBoundaryIntegralDiagnostic=0;    % keep zero (only used for testing)
  



%% Regularisation parameters
CtrlVar.SpeedZero=1e-4;     % needs to be larger than 0 but should also be much smaller than any velocities of interest
CtrlVar.EpsZero=1e-10;      % needs to be larger than 0 but should also be much smaller than any effective strain rates of interest
CtrlVar.etaIntMax=1e10 ;    % max value of effective viscosity
CtrlVar.etaIntMin=1e-6;     % min value of effective viscosity
CtrlVar.Czero=1e-10;        % 
CtrlVar.CAdjointZero=CtrlVar.Czero; % used as a regularisation parameter when calculating dIdCq
CtrlVar.dbdxZero=1;   % when calculating basal shear stresses in the hybrid approximation, a very large bed slope can lead to occillations,
CtrlVar.dbdyZero=1;   % therfore surpress bedslopes to less than 45 degrees

%% Constraints on viscosity and slipperiness
% These constraints are always enforced, but only really of any importance when invertion for A and/or C.
% (Using SIA or the hybrid approximation Cmin MUST be set to 0, or at least to a value much less thatn Czero!)
%
switch lower(CtrlVar.FlowApproximation)
    case 'sstream'
        CtrlVar.Cmin=1e-6;          % a reasonable lower estimate of C is u=C tau^m with min u=1 m/a and max tau=100 dPa => C=u/tau^m=1e-6
    otherwise
        CtrlVar.Cmin=0;          % a reasonable lower estimate of C is u=C tau^m with min u=1 m/a and max tau=100 dPa => C=u/tau^m=1e-6
end
CtrlVar.Cmax=1e10;
CtrlVar.AGlenmin=100*eps;
CtrlVar.AGlenmax=1e10;

%% Non-linear iteration-loop parameters
% The non-linear system is considered solved once the residuals are smaller than NLtol,
% and the normalised chances in u,h and \lambda smaller than du, dh and dl.
%
% The most (arguably even the only) important number is NLtol.
% NLtol is a tolerance on the norm of the actual solution residuals, i.e. the resulting residuals once the solution is
% plugged back into the equation. So NLtol should ALWAYS be set to a small value (for example <1e-10)
%
% The CtrlVar.du/dh/dl are tolerances for the chance in u,h, and \lambda, respectively, between non-linear iterations.
% Although one would expect these to go to zero with increasing iteration number, these are not very reliable
% estimates of the actual error.  Generally set du and dh to not too large value, but do not focus too much on those numbers
% (The error in solving the boundary conditions is always checked internally.)
CtrlVar.NLtol=1e-15; % this is the square of the error, i.e. not root-mean-square error
CtrlVar.du=1e-2;     % tolerance for change in (normalised) speed
CtrlVar.dh=1e-2;     % tolerance for change in (normalised) thickness
CtrlVar.dl=100;      % tolerance for change in (normalised) lambda variables used to enforced BCs
%    Note: there is no need to put any constrains on the Lagrange variables
%    used to enforce the BCs because 1) the BCs are currently always linear,
%    and 2) it is always checked internally that the BCs have been solved correctly.
%    In fact it can be a bad idea to enforce a limit on this chance.
%    (Sometimes the change in lambda between non-linear iteration steps is just a
%    direct response to how the primary variables (u,v,h) change.  The norm
%    of these changes can then be large despite the BCs being exactly fulfilled.)

CtrlVar.NR=1;             % 1 gives Newton-Raphson (use Newton-Raphson whenever possible)
CtrlVar.Piccard=0;        % 1 gives Piccard iteration
CtrlVar.NRviscosity=1;    % if 1 derivatives with respect to viscosity are included in the NR method
CtrlVar.NRbeta2=1;        % if 1 derivatives with respect to slipperiness are included in the NR method
CtrlVar.NRitmax=50;       % maximum number of NR iteration
CtrlVar.Piccarditmax=30;  % maximum number of Piccard iterations
CtrlVar.iarmmax=10;       % maximum number of backtracking steps in NR and Piccard iteration
CtrlVar.NRitmin=1;        % minimum number of NR iteration
CtrlVar.NewtonAcceptRatio=0.5;  % accepted reduction in NR without going into back-stepping
CtrlVar.NewtonBacktrackingBeta=1e-4;  %  affects the Amarijo exit criteria in the back-stepping
CtrlVar.LineSeachAllowedToUseExtrapolation=1; % If true, backtracking algorithm may start with an extrapolation step.
CtrlVar.BacktrackingGammaMin=1e-10;  % smallest step-size in Newton/Piccard backtracking

%% Number of integration points
% if left empty, the number of integration points is set automatically

CtrlVar.niph=[] ;  % number of integration points for uvh in implicit runs, and for the h-solver in semi-implicit runs
CtrlVar.nip=[] ;   % number of integration points for the uv solver


%% Backtracking parameters  -line search 
% Parameters affecting the backtracking algorithm
CtrlVar.BackTrackBeta=0.1 ;               % beta in the Armijo–Goldstein exit condition
CtrlVar.BackTrackMaxIterations=50 ;       % this is plenty
CtrlVar.BackTrackMaxExtrapolations=50  ;  % if set to zero no extrapolation is done (i.e. pure backtracking)
CtrlVar.BackTrackExtrapolationRatio=2.5 ; % ratio between new and old step size in each extrapolation step
CtrlVar.BackTrackMinXfrac=1e-10 ;         % exit backtracking if pos. of minimum is changing by less than this fraction of inital step 
CtrlVar.BackTrackMaxFuncSame=3 ;          % exit backtracking if this many evaluations of cost function resulted in no further decrease of cost funcion
    


%% Lin equation solver parameters


% Linear symmetrical solver is either Matlab \ operator, or Uzawa (outer) iteration
% CtrlVar.SymmSolver can be one of {'Backslash','Uzawa','AugmentedLagrangian'}
CtrlVar.SymmSolver='AugmentedLagrangian';  %

% Linear asymmetrical solver is either Matlab \ operator or Augmented Lagrangian Solver (ALS)
% CtrlVar.AsymmSolver='Backslash';
CtrlVar.AsymmSolver='AugmentedLagrangian';  %
% For asymmetrical indefinite block-structured systems
% the ALS method is almost always better than the default Matlab backslash operator. 
% ALS is an iterative method with an inner and outer iteration. Convergence depends on
% ALSpower. If ALS does not converge then tying a smaller ALSpower
% usually does the trick.
CtrlVar.ALSIterationMin=3;     CtrlVar.ALSIterationMax=25;   CtrlVar.ALSpower=5;  % ALS parameters
CtrlVar.UzawaIterationMin=3;   CtrlVar.UzawaIterationMax=25; CtrlVar.UzawaPower=5;  % Uzawa parameters


CtrlVar.LinSolveTol=1e-10;  % Residual when solving linear system.
                            % If the standard Matlab backslash algorithm is used, default Matlab values apply and this number is not used
                            % For indefinite block-structured systems of the type [A B' ; B 0] [x;y]=[f;g]
                            % the relative residual is defined in standard way as: 
                            % Residual=norm([A B' ; B sparse(m,m)]*[x;y]-[f ; g])/norm([f;g]);   
                            % A value of 1e-10 is arguably an overly small number, in many cases 1e-6 would be considered acceptable





%% Level of information given during a run
% A number of variables affect the information given during a run.
% Generally the higher the number the more information is given.
% 
% Depending on info levels, figures might be plotted as well. However, this is only done
% if corresponding plotting logicals such as CtrlVar.doplots, CtrlVar.doAdaptMeshPlot, etc, are also true.
%
CtrlVar.InfoLevel=1;                     % Overall level of information  
CtrlVar.Report_if_b_less_than_B=0;
CtrlVar.InfoLevelLinSolve=0;
CtrlVar.SymmSolverInfoLevel=0 ;
CtrlVar.InfoLevelAdaptiveMeshing=1;
CtrlVar.InfoLevelAdjoint=100;
CtrlVar.InfoLevelNonLinIt=1;
CtrlVar.InfoLevelBackTrack=1;
CtrlVar.ThicknessConstraintsInfoLevel=1 ;

%  CtrlVar.InfoLevelNonLinIt:
% <0   : no information printed
% >=0  : prints basic convergence information at end of non-linear step
% >=1  : detailed info on residuals given at the end of non-linear step
% >=2  : info on backtracking step as well
% >=10 : calculates/plots additional info on residuals as a function of step size within line search, and rate of convergence
% >=100 : plots residual vectors
CtrlVar.InfoLevelCPU=0;  % if 1 then some info on CPU time usage is given

CtrlVar.StandartOutToLogfile=false ; % if true standard output is directed to a logfile
% name of logfile is  $Experiment.log



%% Adjoint variables
CtrlVar.AdjointGrad='C';  % {'C'|'A'}
CtrlVar.MaxAdjointIterations=1;
CtrlVar.AdjointRestart=0;
CtrlVar.AdjointWriteRestartFile=1;
CtrlVar.NameOfAdjointRestartFiletoWrite='AdjointRestart.mat';
CtrlVar.NameOfAdjointRestartFiletoRead=CtrlVar.NameOfAdjointRestartFiletoWrite;
CtrlVar.NameOfFileForSavingSlipperinessEstimate='C-Estimate.mat';
CtrlVar.NameOfFileForSavingAGlenEstimate='AGlen-Estimate.mat';

CtrlVar.AdjointInitialSearchStepSize=[]; % initial guess for step size in line-search
                                         % If left empty, step size is based on InvStartValues, unless in a restart run when the last converged value is used.
                                         % Mostly usefull for resetting step size in a restart run


% There are a number of different minimisation methods implemented
% currently the QuasiNewtonInversion seems to work best but this may depend on
% the particular case in question.
% 
CtrlVar.AdjointMinimisationMethod='QuasiNewtonInversion';  % works
%CtrlVar.AdjointMinimisationMethod='QuasiNewtonInversion:HessianGuesstimate';  % here an educated guess for 
                                                                              % the Hessian is used
                                                                              % This is not guaranteed to work
                                                                              % because the guesstimate might be
                                                                              % badly wrong. However,
                                                                              % experience has shown that in many cases
                                                                              % this works well and in fact often 
                                                                              % increases the rate of convergence 

%CtrlVar.AdjointMinimisationMethod='FixPointEstimationOfSlipperiness';  % works
%CtrlVar.AdjointMinimisationMethod='AdjointProjectedGradient' ;  % works
%CtrlVar.AdjointMinimisationMethod='MatlabConstrainedMinimisation'; % works
% CtrlVar.AdjointMinimisationMethod='ProjectedBFGS';  % broken


CtrlVar.RescaleAdjointGradient=0;  % rescales analytical gradient to agree with a numerically calculated one (only use for testing purposes)
CtrlVar.CalcBrutForceGradient=0;

CtrlVar.AdjointfScale=1;
CtrlVar.AdjointxScale=1;

CtrlVar.MisfitFunction='uvintegral';
CtrlVar.AdjointGradientEvaluation='integral';
%CtrlVar.MisfitFunction='uvdiscrete';
CtrlVar.NormalizeWithAreas=1 ;  % Cost function normalized with element areas. 
                                % (generally a good idea as it makes gradient indpended of ele size)

CtrlVar.AdjointEleConst=1 ;


CtrlVar.isBarrierC=1     ; CtrlVar.muBarrierCmin=1e-10     ; CtrlVar.muBarrierCmax=1e-10 ;  % note: in most cases with constraints the muBarrier parameters should initially be set to fairly large values
CtrlVar.isBarrierAGlen=1 ; CtrlVar.muBarrierAGlenmin=1e-10 ; CtrlVar.muBarrierAGlenmax=1e-10 ;
CtrlVar.isRegC=1; CtrlVar.isRegAGlen=1;

CtrlVar.RegAGlenMultiplier=1; CtrlVar.RegCMultiplier=1   ;% the regularisation terms are muliplied by these numbers,
% good for increasing/decreasing
% the relative size of the regularisation term

CtrlVar.MisfitMultiplier=1;   % the misfit term is multiplied with this number
                              % (increasing this number makes other terms in the cost funcition (regularisation,barrier)
                              % less important in comparision to the data misfit term.)



% Conjugent gradient parameters: 
CtrlVar.AdjointConjugatedGradients=1; 
CtrlVar.ConjugatedGradientsRestartThreshold=0.2;
CtrlVar.ConjugatedGradientsUpdate='PR'; % (FR|PR|HS|DY)
                                        % FR ;Fletcher-Reeves
                                        % PR :Polak-Ribi\`ere
                                        % HR: Hestenes-Stiefel
                                        % DY :Dai-Yan
                            
CtrlVar.AdjointGradientEleAverage=0;
CtrlVar.AGlenAdjointZero=1e-10; 
CtrlVar.AdjointEpsZero=1e-10;
CtrlVar.AdjointMaxLineSearchIterations=20;
CtrlVar.CalcBruteForceGradient=0;
CtrlVar.RescaleAdjointGradient=0;


CtrlVar.ReadInitialSlipEstimateForInversionFromFile=[];   % If starting inversion with values from a file, give file name for input file with Cest and m, leave empty otherwise
CtrlVar.ReadInitialAGlenEstimateForInversionFromFile=[];  % If starting inversion with values from a file, give file name for input file with AGlenest and n, leave empty otherwise


% BFGS parameters
CtrlVar.Maximum_Number_of_BFGS_updates=250;

CtrlVar.SaveSurfaceData=0;
CtrlVar.UseSyntheticData=0;
CtrlVar.SurfaceDataFile='ForwardC-SlipperyPert.mat'; CtrlVar.SurfaceDataFileLoaded=0;

%CtrlVar.SaveSurfaceData=0; % overridng

CtrlVar.CompareWithAnalyticalSolutions=0;
CtrlVar.CompareResultsWithPreviouslyObtainedResults=0;


%% Numbering of nodes and elements
CtrlVar.sweep=1;              % renumber nodes using a `sweep' plane
CtrlVar.SweepAngle=0.01;      % angle of sweep plane with respect to x axis
CtrlVar.CuthillMcKee=0;       % renumber nodes using sparse reverse Cuthill-McKee ordering

%% Creation of a dumpfile
% Mainly used for testing purposes, but can in principle also be used to generate output data files.
CtrlVar.WriteDumpFile=0;                      % a dumpfile is created containing all variables
CtrlVar.WriteDumpFileStepInterval=1000;       % number of time steps between writing a dump file
CtrlVar.WriteDumpFileTimeInterval=0;          % time interval between writing a dump file

%%  Call a user output routine

CtrlVar.UaOutputsDt=0; % model time interval between calling UaOutputs.m
                       % if set to zero UaOutputs is called at every time step
                       % if set to a negative number, or NaN, UaOutputs is never called
CtrlVar.UaOutputsMaxNrOfCalls=NaN;  % maximum nr of calls to UaOutputs
                                    % Once this limit is reached, the run stops. (Setting this to 1 or some low number
                                    % can sometimes be usefull for testing/control purposes)
                                    % if set to NaN implies no limit to the number of calls 
%% Meshing
% There are various ways of meshing the computational domain.
%
% The simplest option tends to be to use the external mesh generator `gmsh' directly from within Úa
% If that is done, only the outlines of the mesh need to be defined in the variable 'MeshBoundaryCoordinates'
% and the variables CtrlVar.MeshSizeMin and CtrlVar.MeshSizeMax set. However there are various other options.
%
%  For examples of how to generate different type of meshes look at `ExamplesOfMeshGeneration.m'
%
%
% The main Options are:
%
% # Read the initial FE mesh (coordinates, connectivity) directly into Ua from a previously generated .mat file:
%   (set CtrlVar.ReadInitialMesh to true)
% # Use `mesh2d.m' to generate the mesh (set CtrlVar.MeshGenerator='mesh2d')
% # Use `gmesh'  to generate the mesh (set CtrlVar.MeshGenerator='gmesh')
%
% Note: mesh2d is excellent, but unfortunately does not work with matlab2011b or later.
% Therefore `gmesh' is at the moment the best option to use.
%
% In general using `gmesh' involves:
% *             a) creating an inputfile for gmesh (.geo)
% *             b) calling gmesh for that input file
% *             c) reading the resulting gmesh output file (.msh) with the mesh
% All, or some of these three steps can be done withing Úa.
%
% When using gmesh it is possible to select between these option:
%
% *    i)  Directly read existing gmesh output file (.msh)
% *   ii)  First run gmesh with an existing gmesh input file (.geo) and then read the resulting gmesh output file (.msh)
% *   iii) First generate gmesh input file (geo), then run gmesh for that input file, and finally read the resulting gmesh output file (.msh)
%
% To select between 3i and 3iii set CtrlVar.GmeshMeshingMode={'load .msh','mesh domain and load .msh file','create new gmesh .geo input file and mesh domain and load .msh file'}
%
% In most cases running gmesh from within Úa is the simplest way of meshing the domain.
% Only the outlines of the domain need to be defined as well as a few parameters that control
% the element sizes.
%

CtrlVar.ReadInitialMesh=0;    % if true then read FE mesh (coordinates, connectivity) directly from a .mat file 
                              % unless the adaptive meshing option is used, no further meshing is done.
CtrlVar.ReadInitialMeshFileName='ExistingMeshFile.mat';
CtrlVar.SaveInitialMeshFileName='NewMeshFile.mat';
CtrlVar.MeshGenerator='gmesh';  % possible values: {mesh2d|gmesh}

% CtrlVar.GmeshMeshingMode='load .msh'                                                              % option 3i
% CtrlVar.GmeshMeshingMode='mesh domain and load .msh file'                                         % option 3ii
CtrlVar.GmeshMeshingMode='create new gmesh .geo input file and mesh domain and load .msh file';     % option 3iii


CtrlVar.GmeshFile='GmeshFile';  % name of gmesh input/output files (no file extensions)

CtrlVar.GmeshMeshingAlgorithm=1;    % see gmsh manual
                                    % 1=MeshAdapt
                                    % 2=Automatic
                                    % 5=Delaunay
                                    % 6=Frontal
                                    % 7=bamg
                                    % 8=DelQuad (experimental)
                                    

CtrlVar.GmeshBoundaryType='lines';   % (spline|lines)
CtrlVar.GmeshCharacteristicLengthExtendFromBoundary=0;
CtrlVar.GmeshCharacteristicLengthFromCurvature = 0 ;
CtrlVar.GmshGeoFileAdditionalInputLines{1}='   ';  % these lines are added to the gmsh .geo input file each time such a file is created

CtrlVar.OnlyMeshDomainAndThenStop=0; % if true then only meshing is done and no further calculations. Useful for checking if mesh is reasonable
CtrlVar.AdaptMeshAndThenStop=0;

%% Controlling element sizes
% if no adaptive meshing is used then the element size is given by
CtrlVar.MeshSize=10e3;                       % over-all desired element size (however if gmsh is used without automated mesh adapting
                                             % then only CtrlVar.MeshSizeMin and CtrlVar.MeshSizeMax are used)
                                             % 

% if adaptive meshing is used then the range of min to max element sizes is:
CtrlVar.MeshSizeMin=0.1*CtrlVar.MeshSize;    % min element size
CtrlVar.MeshSizeMax=CtrlVar.MeshSize;        % max element size

CtrlVar.MaxNumberOfElements=100e3;           % max number of elements. If #elements larger then CtrlMeshSize/min/max are changed
CtrlVar.MaxNumberOfElementsUpperLimitFactor=1.3;  % if actual number of elements is larger than CtrlVar.MaxNumberOfElements by this factor
                                                  % the domain is remeshed by modifying MeshSizeMin 
CtrlVar.MaxNumberOfElementsLowerLimitFactor=0.0;


%% Pos. thickness constraints,          (-active set-)
% A minimum ice thickness can be enforced in different ways using the following methods:
%  1) `reset method' : simply resetting the thickness to min thickness at node where thickness is less than a prescribed value.
%  2) `active-set' method.
%  3) `thickness-barrier' method
%
% The active-set method is the preferred option and is arguably the only correct way of enforcing min ice thickness.
% The active-set method should therefore be used whenever possible.
% However, in some cases the active set method does not converge, in which case
% options 1) or 3), or combinations thereof, must be used.  If the differences between approach 1) and 2) are small, then using 1)
% allows for shortest computation times
%
% The thickness-barrier method introduces a fictitious surface mass balance term.
% The thickness-barrier method can be used on its own, but should primarily be used in combination with the active-set method
% to improve convergence.
%


CtrlVar.ThickMin=1;                      % minimum allowed thickness without (potentially) doing something about it

% reset method, option 1
CtrlVar.ResetThicknessToMinThickness=0;  % set to 1 to reset thickness values less than ThickMin to ThickMin at each time step (Option 1, not recommended)
CtrlVar.ResetThicknessInNonLinLoop=0;    % if true, thickness in the non-linear iteration of the uvh implicit approach
                                         % is set to zero, provided CtrlVar.ResetThicknessToMinThickness is also true (usually not a good idea)


% active-set method, option 2 
CtrlVar.ThicknessConstraints=1;             % set to 1 to use the active-set method (Option 2, the recommended option).
CtrlVar.ThicknessConstraintsItMax=10  ;     % maximum number of active-set iterations.
                                            % if the maximum number of active-set iterations is reached, a warning is give, but
                                            % the calculation is not stopped. (In many cases there is no need to wait for
                                            % full convergence of the active-set method for each time step.)
                                            % if set to 0, then the active set is updated once and then proceed to next time step.
CtrlVar.ThicknessConstraintsLambdaPosThreshold=0;  % if Thickconstraints are larger than this value they are inactivated, should be zero
CtrlVar.NumberOfActiveThicknessConstraints=0;      % The number of active thickness constraints (just for information, always set initially to zero)
CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints=1000 ; %

% thickness barrier, option 3
CtrlVar.ThicknessBarrier=0;                   % set to 1 for using the barrier method  (Option 3)
CtrlVar.ThicknessBarrierThicknessScale=CtrlVar.ThickMin;     %
CtrlVar.ThicknessBarrierDiagonalFraction=1;   % size of barrier term in comparison to mean abs of diagonal elements
CtrlVar.ThicknessBarrierMinThickMultiplier=2; % exp. barrier is 1 at ThickMin * MinThickMuliplier
CtrlVar.ThicknessBarrierAccumulation=0.01;

%% Advance/Retreat mesh
% This option allows for deactivation/activation of elements based on ice thickness.
% A `background' FE mesh is required. In most cases this background FE mesh will simply be the inital FE mesh
% used at the start of the calculation.
% For advancing glaciers this option must be combined with the active-set method (set CtrlVar.ThicknessConstraints=1)
%
% Note: When combined with the active-set method
% then CtrlVar.ThickMinDeactivateElements must be >= CtrlVar.ThickMin.
% It is usually good for this value to be slightly larger than CtrlVar.ThickMin
% Nodes are only included in the active set if thickness at a node < CtrlVar.ThickMin.
% If CtrlVar.ThickMinDeactivateElements>CtrlVar.ThickMin elements will be eliminated
% before all nodes of that element have (potentially) been included in the active set.
% This reduces somewhat the number of nodes in the active set from what it would otherwise be
% if CtrlVar.ThickMinDeactivateElements=CtrlVar.ThickMin.
%
% Note: Elements are only inactivated if ALL nodes have thickness <CtrlVar.ThickMinDeactivateElements
% Elements are activated once at least one of the nodes leaves the active set. It is therefore
% possible to have a new element where some of the nodes have thickness larger than CtrlVar.ThickMin.
% but less than CtrlVar.ThickMinDeactivateElements.
CtrlVar.FEmeshAdvanceRetreat=0;     % activates the Advance/Retreating mesh option
CtrlVar.FEmeshAdvanceRetreatDT=0.5; % activaton/deactivation done at this time interval
                                    % for CtrlVar.FEmeshAdvanceRetreatDT=0 the activation/deactivation is done
                                    % at every time step (in many cases the best approach)
CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName='BackgroundMeshfile.mat'; % This file is needed for the advance/retreat option
                                                                             % It must contain the variable `MUA_Background'
                                                                             
CtrlVar.ThickMinDeactivateElements=1.01*CtrlVar.ThickMin;% Elements where thickness at all nodes is less than this value are deactivated
CtrlVar.SelectElementsToDeactivateAlgorithm=1; % (1|2)  There are two different methods implemented for selecting
                                                % elements to be deactivated in conjunction with the FEmeshAdvanceRetreat option:
                                                % 1: Eliminate an element if none of the nodes of that element belong to an element
                                                %    where any of the nodal thicknesses are greater than CtrlVar.ThickMinDeactivateElements
                                                % 2: Eliminate elements where all the nodes are less or equal to CtrlVar.ThickMinDeactivateElements
                                                % Method 1 eliminates less (or equal) number of elements than Method 2.
                                                % Method 2 is OK for retreating cases and some advancing cases but can fail                                             
                                                % if the advance is `too quick'.

CtrlVar.MinSurfAccRequiredToReactivateNodes=0;  % If surface accumulation is larger than this, then a node is considered to have positive ice thickness
                                                % and not eliminated.  This is important in cases where, for example, with time the ELA drops down below 
                                                % the top of a mountain peak that is not included in the current FE-mesh.
                                                % This allows for the formation of new isolated glaciated areas.
                                                % Although the default value is zero, it is presumably better to set this to a small positive value.

%% Uniform global mesh refinement
% Mesh can be refined at a start of a run or the start of a restart run by subdividing all triangles into four
% can be useful, for example, for an error estimation
CtrlVar.RefineMeshOnRestart=0;
CtrlVar.RefineMeshOnStart=0;

%% Global adaptive mesh refinement  , adapt mesh 
% There are various adapt meshing options.
% The most general one is global remeshing using explicit error estimate
%
% Global remeshing can be based on one or more of the following
% RefineCriteria:
%           'effective strain rates'
%           '|dhdt|'
%           '||grad(dhdt)||'
%           'dhdt curvature'
%           'thickness gradient'
%           'thickness curvature'
%           'flotation'
%           'f factor'
% the criteria can be combined.
% When two or more criteria are combined RefineCriteria is given as a cell array
%
% The relative importance of different RefineCriteria can be specified by defining `RefineCriteriaWeights'
% These weights affect how small the smallest element will be for a given refinement criteria.
%
% If CtrlVar.RefineCriteriaWeights=1, the whole range CtrlVar.MeshSizeMin to CtrlVar.MeshSizeMax is used
%
% If CtrlVar.RefineCriteriaWeights=0.5 element size will range from
% CtrlVar.MeshSizeMax down to CtrlVar.MeshSizeMin+(CtrlVar.MeshSizeMax-CtrlVar.MeshSizeMin)*(1-RefineCriteriaWeight)
%
% If CtrlVar.RefineCriteriaWeights=0 the criterion is effectively ignored
%
% Examples:
%
%  CtrlVar.RefineCriteria='effective strain rates';  % specifies 'effective strain rates' as the only criterion
%  CtrlVar.RefineCriteriaWeights=[1];                % with a relative weight of unity
%
%  CtrlVar.RefineCriteria={'flotation','||grad(dhdt)||','dhdt curvature','thickness curvature'}; % several criteria used
%  CtrlVar.RefineCriteriaWeights=[1,1,1];  %
%
% In addition the refinement can be limited to an area within a given flotation distance by defining
% CtrlVar.RefineCriteriaFlotationLimit
% For example, for CtrlVar.RefineCriteriaFlotationLimit=[100,NaN]
% the first refinement criteria will only be applied to area where the glacier bed (b) is within 100 vertical distance units from flotation
%

CtrlVar.AdaptMesh=0;          % true if adapt meshing is used, no remeshing is done unless this variable is true
CtrlVar.MeshRefinementMethod='explicit:global';    % can have any of these values:
                                                   % 'explicit:global'
                                                   % 'explicit:local'
                                                   % 'implicit:global'  (broken at the moment, do not use)
                                                   % 'implicit:local'   (broken at the moment, do not use)

                                                   
CtrlVar.AdaptMeshInitial=1  ; % remesh in first iteration (Itime=1)  even if mod(Itime,CtrlVar.AdaptMeshInterval)~=0.
CtrlVar.AdaptMeshInterval=1 ; % remesh whenever mod(Itime,CtrlVar.AdaptMeshInterval)==0

CtrlVar.AdaptMeshAndThenStop=0;    % if true, then mesh will be adapted but no further calculations performed
                                   % useful, for example, when trying out different remeshing options (then use CtrlVar.doAdaptMeshPlots=1 to get plots)


                                           
CtrlVar.hpower=1;         % used to go from an error estimate to a size estimate for an element
                          % h=1/error^hpower ,  where `h' is the desired element size and `error' a local
                          % error estimate

CtrlVar.AdaptMeshIterations=1;  % Number of global adapt mesh iterations within each adapt step
                                % This is seldom anything else but 1, except potentially in a either diagnostic calculation (time independent)
                                % or at the start of a transient (prognostic) calculation where the initial mesh is very coarse
                                % Note that when using gmsh the mesh refinement is always based on indicators given at the nodal points of the original mesh
                                % therefore using a few AdaptMeshIterations can sometimes be be a good idea.

                                
CtrlVar.LocalAdaptMeshSmoothingIterations=5;  % Number of Laplace mesh smoothing iterations used in local mesh refinement
CtrlVar.LocalAdaptMeshRatio=0.25;             % The maximum number of elements subdivided during each local mesh refinement step
                                              % as a fraction of the total number of elements.

CtrlVar.MaxRatioOfChangeInEleSizeDuringAdaptMeshing=5;   % put a strickt limit on how much ele sizes change during single
CtrlVar.MinRatioOfChangeInEleSizeDuringAdaptMeshing=1/5; % adaptive meshing step to avoid excessive changes.
                                                         % This does not apply to local mesh refinement where in each adapt step
                                                         % the elements are always only refined, and then always by a factor of two.

CtrlVar.RefineCriteria={'flotation','||grad(dhdt)||','dhdt curvature'};
CtrlVar.RefineCriteriaWeights=[0.1,1,1];                %  
CtrlVar.RefineCriteriaFlotationLimit=[NaN,NaN,NaN];     % Refine criteria is only applied to elements which are this close to flotation
                                                        %        useful to restrict refinement to an area in the (vertical) vicinity of the grounding line
                                                        %        Set to NaN if to be ignored and applied to all regions irrespective of how close to flotation

CtrlVar.NumberOfSmoothingErrorIndicatorIterations=1;    % each of the error indicators can be smooth over neighbouring elements
                                                        % this is done by calculating average values for each element based its nodal values
                                                        % and then interpolating from elements to nodes using number of elements that a node is
                                                        % attached to as a weighting factor.  This introduces a smoothing that is related to
                                                        % connectivity as opposed to spatial distance.
                                                        % This kind of smoothing is never done for the 'flotation' and the `f factor' cases
                                                        % as the spread/smoothing can be determined directly by CtrlVar.RefineDiracDeltaWidth

CtrlVar.RefineDiracDeltaWidth=100;  % for `flotation' and 'f factor' the zone within this vertical distance from flotation is refined
CtrlVar.RefineDiracDeltaOffset=0;   %
                              
                                

%% `Time geometries' are boundary geometries that change with time.
% This can be used, for example, to simulate a calving event.
CtrlVar.TimeGeometries.Flag=0;             % true if domain geometry is changed during the run (e.g prescribed calving event)

%% Mesh morphing:
%
% (mesh morphing around a moving grounding line is currently broken. This 
% looked like a good idea, but really is only going to work if the grounding line has
% a simple shape and the topology of the grounding lines does not change.
% Basically a too limited option for practical use.)
CtrlVar.MeshMorphing=0;       % true for mesh-morphing where the mesh is morphed onto moving grounding line
                              % this is a very elegant method, but only works if there is just one grounding line
                              % and therefore not really that useful in a general 2HD situation.
CtrlVar.GLmeshing=0;          % GL meshing based on morphing
CtrlVar.GLtension=1;          % tension of spline used in GL morphing, 1: no smoothing; 0: straight line
CtrlVar.GLds=CtrlVar.MeshSizeMin ; % edge length along GL when using GL meshing




%% Parameters affecting floating/grounded mask

CtrlVar.kH=1;   % kH -> infty gives an exact Heaviside and delta functions.
                % kH=1 implies a grounding line "width" of 1 m up and down from floating condition
                % kH=10 implies a grounding line "width" of 1/10 m up and down from floating condition
CtrlVar.Hh0=0;  % offset is Heaviside function when calculating GF field

%% Parameters affecting calculation of grounding line
% The grounding line position does not enter any calculations done by Úa. 
% The grounding line is primarirly calculated for plotting purposes.
CtrlVar.GLthreshold=0.5;  % used to define position of GL with respect to the values of the Heaviside function (1 fully grounded, 0 fully floating)
CtrlVar.GLsubdivide=0;    % If 0/false the grounding line is determined based on GL.node values at corners only (using GLthreshold). If 1/true
                          % then all nodal values of 6-node and 10-node triangles are also used. This is done by splitting those into 4 and 9 triangles, respectivily


%%

CtrlVar.CalvingFrontFullyFloating=0;  % if true then the natural BC is only covers a freely floating calving front (do not change, only for testing)
CtrlVar.GroupRepresentation=0;

CtrlVar.AGlenisElementBased=0;
CtrlVar.CisElementBased=0;

CtrlVar.registerGLpos=0;
CtrlVar.Parallel=0; CtrlVar.ParallelAssembly=0;
CtrlVar.registerGLposFilename='GLposition';
CtrlVar.calcDerivedGLquantitiesForGLeleOnly=1;


%% Adaptive Time Stepping Algorithm (ATSA)   (adapt time step)
% The adaptive-time-stepping algorithm is based on the idea of keeping the number of non-linear iterations
% close to a certain target (CtrlVar.ATSTargetIterations).
% This is a simple but highly effective method.  However, as the ATSA is not based on any kind of error estimates,
% it does not guarantee that errors will not accumulate, etc, and ATSA does not work for linear problems.
%
% The main idea is to aim at a time step that limits the number of non-linear iteration to a relatively small number.
% If so, then most likely the Newton-Raphson iteration is in the quadratic regime.
% Experience has shown that a target number of iterations (CtrlVar.ATSTargetIterations) within 3 to 5 is good for this purpose
%
% Time step is increased if r<1 where
%
%     r=N/M
%
%   where
%   N is the max number of non-linear iteration over last n time steps
%   M is the target number of iterations
%
%   here
%     M=CtrlVar.ATSTargetIterations
%   and
%     n=CtrlVar.ATSintervalUp
%
%   (N does not need to be specified.)
%
%  Currently the time step is only decreased if either:
%        a) the number of non-linear iterations in last time step was larger than 25
%        b) number of iterations over last n times steps were all larger than 10
%  where n=CtrlVar.ATSintervalDown
%
% There are some further modifications possible:
%  -time step is adjusted so that time interval for making transient plots (CtrlVar.TransientPlotDt) is not skipped over
%  -time step is not increased further than the target time step CtrlVar.ATStimeStepTarget
%  -time step is adjusted so that total simulation time does not exceed CtrlVar.TotalTime
%
%
%
CtrlVar.AdaptiveTimeStepping=1 ;    % true if time step should potentially be modified
CtrlVar.ATStimeStepTarget=1000.0;   % maximum time step size allowed
CtrlVar.ATStimeStepFactorUp=2 ;     % when time step is increased, it is increased by this factor
CtrlVar.ATStimeStepFactorDown=10 ;  % when time step is decreased, it is decreased by this factor
CtrlVar.ATSintervalUp=5 ;           %
CtrlVar.ATSintervalDown=3 ;         %
CtrlVar.ATSTargetIterations=4;      % if number of non-lin iterations has been less than ATSTargetIterations for
                                    % each and everyone of the last ATSintervalUp iterations, the time step is
                                    % increased by the factor ATStimeStepFactorUp
                                    
                                    
%% Mass-balance geometry feedback
% If the mass balance is a function of geometry, an additional non-linearity is introduced to transient runs.
% This non-linearity can be solved in a fully consisten way using the Newton-Raphson method provided the user
% supplies the gradient of the mass balance with respect to thickness.
%
CtrlVar.MassBalanceGeometryFeedback=0;  % If the mass balance depends on geometry then
                                        % setting this parameter to either 1, 2 or 3 has the effect of 
                                        % the mass-balance beeing updated within the non-linear transient-loop.
                                        % In principle this parameter should always be set to 3, but in practice the
                                        % feedback can often be sufficiently well acounted for by simply updating 
                                        % mass balance at each and every time step (i.e. option 0).
                                        %  
                                        %  0 : no mass-balance geometry feedback considered within non-lin iteratin loop 
                                        %      (however, as always, mass balanced is updated at each time step)
                                        %  1 : mass-balance feedback included at the start of each non-lin iteration, 
                                        %      but not within the backtracking step.
                                        %  2 : Feeback included in non-lin loop, both at the beginning of each NR iteration, 
                                        %      and within backtracking step.  
                                        %  3 : Consistent mass-balance feedback algorithim. As option 2, but with 
                                        %      the gradient of the mass-balance with respect to thickness added to the
                                        %      left-hand side of the NR system. Requires the user to supply this gradient through 
                                        %      `DefineMassBalance.m'. Doing so can lead to a drastic reduction 
                                        %      in the number of NR steps required.
                                        %
                                        %  If mass balance depends on geometry then always try to use:
                                        %  CtrlVar.MassBalanceGeometryFeedback=3
                                        %  and then also give dasdh and dabdh as return arguments in DefineMassBalance.
                                        %
                                        %  However, if there is no dependency of the mass balance on geometry, always set
                                        %  CtrlVar.MassBalanceGeometryFeedback=0
                                        %  as doing so avoids calls to DefineMassBalance.m within the non-line loop.
                                        %
CtrlVar.MassBalanceGeometryFeedbackDamping=0;  % Dampens the update in surface mass balance.
                                               % If not equal to zero, then the actual mass-balance value used at the end of the time step,
                                               % becomes a weighted average of that at the beginning and the (correct) value at the 
                                               % end of the time step.
                                               % The value must be in the range [0,1]
                                               % Only use this if encountering convergence problems.  
                                               % Should always be equal to 0 if only possible.
                                               % If not equal to 0, the algorithim converges to a wrong solution (!),
                                               % although the error might be very small if mass-balance geometry feedback is not that strong.
      
%%
CtrlVar.MeltNodesDefinition='Edge-Wise';

%% Mapping from Mesh1 to Mesh2
% when a new FE is created, the values from the old mesh need to be mapped onto the new mesh
% if the boundary of the mesh has not changed this only involves interpolation
% but if the boundary of Mesh2 is different from that of Mesh1 then extrapolation might be involved
% In this case one must check for such points and give them some sensible outside value.
% 
CtrlVar.Mesh1To2CheckForPointsOutsideMesh1AndInsideConvexHull=1 ; % for non evolving mesh boundaries, can be set to 0/false
CtrlVar.InpolyTol=0.1;       % tolerance when checking inside outpoints using the `inpoly' m-file, should be small compared to size of any element

%% internal variables  (do not change these)
CtrlVar.FE2dTransientPlotsCounter=0; 
CtrlVar.MeshChanged=0;               

                           
end



