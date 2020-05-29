
function [RunInfo,varargout]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin)

nVar=length(varargin) ;
varargout=cell(nVar,1);



switch CtrlVar.MapOldToNew.method
    
    
    case "scatteredInterpolant"
        
        
        [RunInfo,varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin{:}) ;
        
    case "FE form functions"
        
        [RunInfo,varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions(CtrlVar,RunInfo,MUAold,MUAnew,varargin{:});
end


%% if testing

if CtrlVar.MapOldToNew.Test
    
    tMapOld=tic;
    
    nVar=length(varargin) ;
    varargout=cell(nVar,1);
    [varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant(CtrlVar,MUAold,MUAnew,OutsideValues,varargin{:}) ;
    tMapOld=toc(tMapOld);
    
    Test=varargout ;
    
    tMapNew=tic;
    varargout=cell(nVar,1);
    [varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions(CtrlVar,MUAold,MUAnew,varargin{:});
    tMapNew=toc(tMapNew);
    
    fprintf(' tMapOld \t \t tMapNew \n %f \t \t %f \n ',tMapOld,tMapNew)
    
    for I=1:nVar
        [norm(varargout{I}-Test{I})  norm(Test{I})]
    end
    
end


end
