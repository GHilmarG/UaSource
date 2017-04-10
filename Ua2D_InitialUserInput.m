function [UserVar,CtrlVar,MeshBoundaryCoordinates]=Ua2D_InitialUserInput(UserVar,CtrlVar)


%%  This input file is used if Ua is run directly from the source code folder.
%  
% You will always need to create your own version of this file and put it into
% you working directory together with the other user-input m-files. 
%
%


%%
UserVar.MisExperiment='ice0';            % This I use in DefineMassBalance
UserVar.Outputsdirectory='ResultsFiles'; % This I use in UaOutputs

%%

CtrlVar.Experiment=['MismipPlus-',UserVar.MisExperiment];   
%% Types of run
%
CtrlVar.TimeDependentRun=1; 
CtrlVar.TotalNumberOfForwardRunSteps=3;
CtrlVar.TotalTime=100;
CtrlVar.Restart=0;  


CtrlVar.dt=0.01; 
CtrlVar.time=0; 

CtrlVar.UaOutputsDt=0; % interval between calling UaOutputs. 0 implies call it at each and every run step.
                       % setting CtrlVar.UaOutputsDt=1; causes UaOutputs to be called every 1 years.
                       % This is a more reasonable value once all looks OK.

CtrlVar.ATStimeStepTarget=1;
CtrlVar.WriteRestartFile=1;

%% Reading in mesh
CtrlVar.ReadInitialMesh=0;    % if true then read FE mesh (i.e the MUA variable) directly from a .mat file
                              % unless the adaptive meshing option is used, no further meshing is done.
CtrlVar.ReadInitialMeshFileName='AdaptMesh.mat';
CtrlVar.SaveInitialMeshFileName='NewMeshFile.mat';
%% Plotting options
CtrlVar.doplots=1;
CtrlVar.PlotMesh=1; 
CtrlVar.PlotBCs=1;
CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=1;
CtrlVar.doRemeshPlots=1;
CtrlVar.PlotXYscale=1000; 
%%

CtrlVar.TriNodes=3;


CtrlVar.NameOfRestartFiletoWrite=['Restart',CtrlVar.Experiment,'.mat'];
CtrlVar.NameOfRestartFiletoRead=CtrlVar.NameOfRestartFiletoWrite;




%% adapt mesh
%CtrlVar.InfoLevelAdaptiveMeshing=100;


% very coarse mesh resolution
CtrlVar.MeshSize=20e3;       % over-all desired element size
CtrlVar.MeshSizeMax=20e3;    % max element size
CtrlVar.MeshSizeMin=0.05*CtrlVar.MeshSize;     % min element size

% reasonably fine mesh resolution
%CtrlVar.MeshSize=8e3;       % over-all desired element size
%CtrlVar.MeshSizeMax=8e3;    % max element size
%CtrlVar.MeshSizeMin=200;    % min element size

CtrlVar.MaxNumberOfElements=250e3;           % max number of elements. If #elements larger then CtrlMeshSize/min/max are changed

CtrlVar.AdaptMesh=1;           % 
CtrlVar.SaveAdaptMeshFileName='AdaptMesh.mat'; 



CtrlVar.AdaptMeshInitial=1 ;       % if true, then a remeshing will always be performed at the inital step
CtrlVar.AdaptMeshAndThenStop=0;    % if true, then mesh will be adapted but no further calculations performed
                                   % useful, for example, when trying out different remeshing options (then use CtrlVar.doRemeshPlots=1 to get plots)
CtrlVar.doAdaptMeshPlots=0;       % if true and if CtrlVar.doplots true also, then do some extra plotting related to adapt meshing

%CtrlVar.RefineCriteria={'flotation','thickness curvature','||grad(dhdt)||'};
%CtrlVar.RefineCriteriaWeights=[1,1,1];                %  
CtrlVar.RefineCriteriaFlotationLimit=[NaN,NaN];     

CtrlVar.RefineCriteria={'flotation','thickness gradient'};
CtrlVar.RefineCriteriaWeights=[1,0.75];                %  

CtrlVar.RefineCriteria={'thickness gradient'};
CtrlVar.RefineCriteriaWeights=[1];                %  
  
CtrlVar.AdaptMeshInterval=1;  % number of run-steps between mesh adaptation
CtrlVar.AdaptMeshIterations=1;


CtrlVar.MeshAdapt.GLrange=[10000 5000 ; 3000 2000];



%% Pos. thickness constraints
CtrlVar.ThickMin=1; % minimum allowed thickness without (potentially) doing something about it
CtrlVar.ResetThicknessToMinThickness=0;  % if true, thickness values less than ThickMin will be set to ThickMin
CtrlVar.ThicknessConstraints=1  ;        % if true, min thickness is enforced using active set method
CtrlVar.ThicknessConstraintsItMax=5  ; 

%%

xd=640e3; xu=0e3 ; yr=0 ; yl=80e3 ;  
MeshBoundaryCoordinates=[xu yr ; xu yl ; xd yl ; xd yr];

 
end
