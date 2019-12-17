
function varargout=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,OutsideValues,varargin)

nVar=length(varargin) ;
varargout=cell(nVar,1);

switch CtrlVar.MapOldToNew.method
    
    
    case "scatteredInterpolant"
        
        
        [varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant(CtrlVar,MUA1,x2,y2,OutsideValues,varargin{:}) ;
        
    case "FE form functions"
        
        [varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions(CtrlVar,MUA1,x2,y2,varargin{:});
end


%% if testing

if CtrlVar.MapOldToNew.Test
    
    tMapOld=tic;
    
    nVar=length(varargin) ;
    varargout=cell(nVar,1);
    [varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant(CtrlVar,MUA1,x2,y2,OutsideValues,varargin{:}) ;
    tMapOld=toc(tMapOld);
    
    Test=varargout ;
    
    tMapNew=tic;
    varargout=cell(nVar,1);
    [varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions(CtrlVar,MUA1,x2,y2,varargin{:});
    tMapNew=toc(tMapNew);
    
    fprintf(' tMapOld \t \t tMapNew \n %f \t \t %f \n ',tMapOld,tMapNew)
    
    for I=1:nVar
        [norm(varargout{I}-Test{I})  norm(Test{I})]
    end
    
end


end
