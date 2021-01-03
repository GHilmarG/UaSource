
function [UserVar,CtrlVar]=GetInitialInputs(UserVar,CtrlVar,varargin)
    
    
    OldInputFile="Ua2D_InitialUserInput.m" ; 
    NewInputFile="DefineInitialInputs.m" ; 
    
    if exist(OldInputFile,'file')==2  && ~(exist(NewInputFile,'file')==2)
        
        warning("OldInputFormat:Ua2D_InitialUserInput","Ua2D_InitialUserInput.m  no longer used. Rename that file to DefineInitialInputs.m")
        
    end
    
    
    InputFile="DefineInitialInputs.m" ;
    TestIfInputFileInWorkingDirectory(InputFile) ;
    
    
    
    % Get user-defined parameter values
    %  CtrlVar,UsrVar,Info,UaOuts
    if nargin("DefineInitialInputs.m")>2
        [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar,varargin{:});
    else
        [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar);
    end
    
    CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
    
    
end
