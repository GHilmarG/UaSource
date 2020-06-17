
function [UserVar,CtrlVar]=GetInitialInputs(UserVar,CtrlVar,varargin)
    
    
    OldInputFile="Ua2D_InitialUserInput.m" ; 
    NewINputFile="DefineInitialInputs.m" ; 
    
    if exist(OldInputFile,'file')==2  && ~(exist(NewINputFile,'file')==2)
        
        warning("OldInputFormat:Ua2D_InitialUserInput","Ua2D_InitialUserInput.m  no longer used. Rename that file to DefineInitialInputs.m")
        
    end
    
    
    
    if ~exist(fullfile(cd,'DefineInitialInputs.m'),'file')
        
        fprintf('The input-file DefineInitialInputs.m not found in the working directory (%s).\n',pwd)
        fprintf('This input-file is required for Ua to run.\n')
        error('Ua2D:InputFileNotFound','DefineInitialInputs.m not found in working directory.')
        
    end
    
    
    
    
    % Get user-defined parameter values
    %  CtrlVar,UsrVar,Info,UaOuts
    if nargin("DefineInitialInputs.m")>2
        [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar,varargin{:});
    else
        [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar);
    end
    
    CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
    
    
end
