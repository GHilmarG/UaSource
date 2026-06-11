
function [RunInfo,varargout]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin)

nVar=length(varargin) ;
varargout=cell(nVar,1);



switch CtrlVar.MapOldToNew.method
    
    
    case "scatteredInterpolant"
        
        
        [RunInfo,varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin{:}) ;
        
    case "FE form functions"
        
        [RunInfo,varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions(CtrlVar,RunInfo,MUAold,MUAnew,varargin{:});

    case "ShapeAndScattered"  

        [RunInfo,varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingShapeAndScattered(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin{:});
        

    otherwise

        error(" case not found")
end




end
