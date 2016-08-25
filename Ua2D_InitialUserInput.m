
function [UserVar,CtrlVar,time,dt,MeshBoundaryCoordinates]=Ua2D_InitialUserInput(CtrlVar)
%%
% [UserVar,CtrlVar,time,dt,MeshBoundaryCoordinates]=Ua2D_InitialUserInput(CtrlVar)
% 
% Set the values of various fields of the run control variable CtrlVar as well
% as the time, time step, and the MeshBoundaryCoordinates
%
%
%
% See also Ua2D_DefaultParameters.m
%
%
%%
    warning('Ua:DefaultDefine','Using default Ua2D_InitialUserInput')
    warning('off','Ua:DefaultDefine')
    UserVar.Experiment='DefaultRun';
   
    CtrlVar.doPrognostic=0 ;
    CtrlVar.doDiagnostic=1  ;
    CtrlVar.doInverseStep=0;
    CtrlVar.Restart=0;  
    
    CtrlVar.InitialDiagnosticStep=1; CtrlVar.Implicituvh=1;
    time=0 ; dt=1; CtrlVar.TotalNumberOfForwardRunSteps=1;
    
    CtrlVar.FlowApproximation='hybrid';
    
    %%
    xd=200e3; xu=-200e3 ; yl=200e3 ; yr=-200e3;
    MeshBoundaryCoordinates=flipud([xu yr ; xd yr ; xd yl ; xu yl]);
    CtrlVar.GmshGeoFileAdditionalInputLines{1}='Periodic Line {1,2} = {3,4};';  % these lines are added to the gmsh .geo input file each time such a file is created
    CtrlVar.OnlyMeshDomainAndThenStop=0;
    %CtrlVar.GmshMeshingAlgorithm=8;  % see gmsh manual 
    CtrlVar.TriNodes=6;   % [3,6,10]
    CtrlVar.MeshSize=200e3;
    CtrlVar.MeshSizeMin=0.001*CtrlVar.MeshSize;
    CtrlVar.MeshSizeMax=CtrlVar.MeshSize;
    
    
    CtrlVar.AdaptMesh=1;
    CtrlVar.AdaptMeshInitial=1  ; 
    CtrlVar.AdaptMeshIterations=4;
    CtrlVar.NumberOfSmoothingErrorIndicatorIterations=1;
    CtrlVar.AdaptMeshAndThenStop=0;
    %CtrlVar.MeshRefinementMethod='implicit'; 
    CtrlVar.MeshRefinementMethod='explicit:local'; 
    CtrlVar.RefineCriteria='effective strain rates';
    %CtrlVar.RefineCriteria='||grad(dhdt)||';
    %CtrlVar.RefineCriteria='dhdt curvature';
    %CtrlVar.RefineCriteria='|dhdt|';
    %CtrlVar.RefineCriteria='thickness curvature';
    %CtrlVar.RefineCriteria='thickness gradient';
    CtrlVar.MaxNumberOfElements=2500;
    
   
    
    %%
    CtrlVar.LineSeachAllowedToUseExtrapolation=1;
    
    %%

    CtrlVar.CompareWithAnalyticalSolutions=0;
    CtrlVar.CompareResultsWithPreviouslyObtainedResults=0;
    
    CtrlVar.WriteRestartFile=1;
    
    CtrlVar.doplots=1;
    CtrlVar.doRemeshPlots=1;
    CtrlVar.doAdaptMeshPlots=1; % if true and if CtrlVar.doplots true also, then do some extra plotting related to adapt meshing
    CtrlVar.PlotNodes=1;       % If true then nodes are plotted when FE mesh is shown
    CtrlVar.PlotXYscale=1000;     % used to scale x and y axis of some of the figures, only used for plotting purposes
    CtrlVar.PlotLabels=0 ; CtrlVar.PlotMesh=1; CtrlVar.PlotBCs=1;CtrlVar.PlotNodes=1;
    CtrlVar.TransientPlotDt=NaN;   % model time interval between calls to `FE2dTransientPlots.m'
    CtrlVar.InfoLevelNonLinIt=1;
    
    
end
