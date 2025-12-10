function [points,edge,part]=MeshBoundaryCoordinates2Mesh2dFormat(CtrlVar,MeshBoundaryCoordinates)


nNaN=numel(find(isnan(MeshBoundaryCoordinates(:,2))));


if CtrlVar.Mesh2dInputFormat==1
    
    if nNaN==0
        
        
        points=MeshBoundaryCoordinates;
        edge=[];
        part=[];
        
    else
        
        
        if ~isnan(MeshBoundaryCoordinates(1,2))
            MeshBoundaryCoordinates=[1 NaN ; MeshBoundaryCoordinates];
        end
        
        if ~isnan(MeshBoundaryCoordinates(end,2))
            MeshBoundaryCoordinates=[MeshBoundaryCoordinates ; NaN NaN];
        end
        
        Inan=isnan(MeshBoundaryCoordinates(:,2));
        
        points=MeshBoundaryCoordinates(~Inan,:);
        edge=points*0;
        iNaN=find(Inan);
        
        for k=1:numel(iNaN)-1
            N1=iNaN(k)-k+1;
            N2=iNaN(k+1)-k-1;
            edge(N1:N2,1)=(N1:N2)';
            edge(N1:N2,2)=[(N1+1:N2)';N1];
        end
        
        part=[];
        
    end
    
else
    
    points=MeshBoundaryCoordinates;
    edge=CtrlVar.Mesh2d.edge;
    part=CtrlVar.Mesh2d.part;
    
end

end