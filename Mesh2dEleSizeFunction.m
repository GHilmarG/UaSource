function Mesh2dEleSize=Mesh2dEleSizeFunction(pp,CtrlVar,UserVar,EleSizeScalarField,F)

persistent FEleSize




if isempty(EleSizeScalarField)
    
    Mesh2dEleSize=zeros(size(pp,1),1)+CtrlVar.MeshSize;
    
else
    
    
    
    x=EleSizeScalarField.xy(:,1);
    y=EleSizeScalarField.xy(:,2);
    v=EleSizeScalarField.EleSize;
    
    if isempty(FEleSize) ||  ~isequal(v,FEleSize.Values)
        
        FEleSize = scatteredInterpolant(x,y,v);

    end

    Mesh2dEleSize=FEleSize(pp(:,1),pp(:,2));
    
    
end




end