function varargout=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,OutsideValues,varargin)

varargout=cell(length(varargin),1);

if CtrlVar.TestMapOldNew
    tMapOld=tic;
    
    [varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant(CtrlVar,MUA1,x2,y2,OutsideValues,varargin{:}) ;
    tMapOld=toc(tMapOld);

    Test=varargout{1} ;
    
    tMapNew=tic;
    [varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions(CtrlVar,MUA1,x2,y2,varargin{:});
    tMapNew=toc(tMapNew);
    
    fprintf(' tMapOld \t \t tMapNew \n %f \t \t %f \n ',tMapOld,tMapNew)
    
    [norm(varargout{1}-Test)  norm(Test)]
    
end


end
