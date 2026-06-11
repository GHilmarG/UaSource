function Mesh2dEleSize=Mesh2dEleSizeFunction(pp,CtrlVar,UserVar,EleSizeScalarField,F)
    
    %%
    %
    % Mesh2dEleSize=Mesh2dEleSizeFunction(pp,CtrlVar,UserVar,EleSizeScalarField,F)
    %
    %   A mesh-size function for the mesh-generator mesh2d
    %
    %   This function is called internally by the mesh-generator mesh2d. 
    %   Consult mesh2d documentation for further details. 
    %    
    %   Note:   The mesh2d documentation (see for example refine2.m) refers
    %           to this function has hfun.m
    %
    %   Mesh2dEleSize : desired element size at x, y locations p(:,1) and p(:,2)
    %
    %   F   :   Ua fields, this will be either be empty, or based on
    %           a previous mesh. 
    %
    %%
    
    persistent FEleSizeInterpolant
    
    xMesh2d=pp(:,1);
    yMesh2d=pp(:,2);
    
    if isempty(EleSizeScalarField)
        
        Mesh2dEleSize=zeros(size(xMesh2d))+CtrlVar.MeshSize;
        
    else
        
        
        
        x=EleSizeScalarField.xy(:,1);
        y=EleSizeScalarField.xy(:,2);
        v=EleSizeScalarField.EleSize;
        
        if isempty(FEleSizeInterpolant) ||  ~isequal(v,FEleSizeInterpolant.Values)
            
            FEleSizeInterpolant = scatteredInterpolant(x,y,v);
            
        end
        
        Mesh2dEleSize=FEleSizeInterpolant(xMesh2d,yMesh2d);
        
        
    end
    
   Mesh2dEleSize=max(CtrlVar.MeshSizeMin,Mesh2dEleSize);
   Mesh2dEleSize=min(CtrlVar.MeshSizeMax,Mesh2dEleSize);
    
    
    
end