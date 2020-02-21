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
clear DefineGeometry
clear DefineInputsForInverseRun
clear DefineDensities
clear DefineDesiredEleSize
clear DefineBoundaryConditions
clear DefineMassBalance
clear UaOutputs
clear GetDesiredEleSize
% Ua utilities
clear PlotMeshScalarVariable
clear PlotFEmesh
clear FindOrCreateFigure

end
