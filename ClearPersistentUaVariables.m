function ClearPersistentUaVariables()


%% Clear any persistent variables
clear AdaptiveTimeStepping
clear AdaptMesh
clear BCs2MLC
clear CostFunctionValueAndGradient
clear EleAverageInterpolate
clear JGH
clear NrOfIntegrationPoints
clear LocalMeshRefinement
clear MeshAdvanceRetreat
clear Mesh2dEleSizeFunction
clear multiWaitbar
clear NewConjugatedGrad
% User input files
clear DefineSlipperyDistribution
clear DefineAGlenDistribution
clear DefineGeometry     % I'm not clearing DefineGeometryAndDensities, maybe good so that interpolants do not need to be reloaded
clear DefineInputsForInverseRun
clear DefineDensities
clear DefineDesiredEleSize
clear DefineBoundaryConditions
clear DefineMassBalance
clear DefineOutputs
clear GetDesiredEleSize
clear DefineCalving
clear DefineCalvingAtIntegrationPoints
% Ua utilities
clear PlotMeshScalarVariable
clear PlotFEmesh
clear FindOrCreateFigure
clear LevelSetEquation

end
