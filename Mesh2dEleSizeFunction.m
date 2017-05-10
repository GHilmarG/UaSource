function HH=Mesh2dEleSizeFunction(pp,CtrlVar,GmshBackgroundScalarField)

persistent F

if isempty(GmshBackgroundScalarField)
    
    HH=zeros(size(pp,1),1)+CtrlVar.MeshSize;
    
else
    
    
    
    x=GmshBackgroundScalarField.xy(:,1);
    y=GmshBackgroundScalarField.xy(:,2);
    v=GmshBackgroundScalarField.EleSize;
    
    if isempty(F) ||  ~isequal(v,F.Values)
        
        F = scatteredInterpolant(x,y,v);

    end

    HH=F(pp(:,1),pp(:,2));
    
    
end




end