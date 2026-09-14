
function [RunInfo,varargout]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin)

nVar=length(varargin) ;
varargout=cell(nVar,1);



switch CtrlVar.MapOldToNew.method
    
    
    case "scatteredInterpolant"
        
        
        [RunInfo,varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin{:}) ;
        
    % case "FE form functions" This had over time gone out of use and is now both broken and redundant
    % 
    %     [RunInfo,varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions(CtrlVar,RunInfo,MUAold,MUAnew,varargin{:});

    case "ShapeAndScattered"  

        [RunInfo,varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingShapeAndScattered(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin{:});
        

    otherwise

        ValidMethods=["scatteredInterpolant","ShapeAndScattered"] ;

        Method=string(CtrlVar.MapOldToNew.method) ;

        if Method=="FE form functions"

            error("MapNodalVariablesFromMesh1ToMesh2:ObsoleteMethod", ...
                "CtrlVar.MapOldToNew.method='%s' is no longer supported and has been removed. \n"+ ...
                "\n"+ ...
                "Use 'ShapeAndScattered' instead. It uses the FE form functions for all nodes \n"+ ...
                "inside the old mesh, and OutsideValues or scattered extrapolation for those \n"+ ...
                "outside, and is therefore a strict superset of the option you have selected. \n"+ ...
                "\n"+ ...
                "Note that MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions still exists. \n"+ ...
                "It is now called directly, as (CtrlVar,MUA,x,y,...), to map nodal values onto \n"+ ...
                "arbitrary (x,y) locations rather than onto the nodes of another mesh. \n", ...
                Method)

        else

            error("MapNodalVariablesFromMesh1ToMesh2:UnknownMethod", ...
                "CtrlVar.MapOldToNew.method='%s' is not a recognised option. \n"+ ...
                "Valid options are: %s \n", ...
                Method,join("'"+ValidMethods+"'",", "))

        end

end




end
