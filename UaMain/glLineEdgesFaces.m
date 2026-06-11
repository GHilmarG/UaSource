function [MeshBoundaryCooWithGLcoo,edge,face,xGL,yGL]=glLineEdgesFaces(GF,coordinates,connectivity,MeshBoundaryCoordinates,CtrlVar)
    
    % returns edges and faces to be used in `meshfaces' based on position of grounding line
    % currently only works one GL within domain (can be fairly easily modified though)
    %
    
    [~,~,xGL,yGL]=GLgeometry(connectivity,coordinates,GF,CtrlVar);
    
    if any(isnan(xGL))
        ind=find(isnan(xGL))-1;
        xGL=xGL(1:ind(1)); yGL=yGL(1:ind(1));
    end
    % equally spaced and possible smoothed
    
    [xGL,yGL] = Smooth2dPos(xGL,yGL,CtrlVar);

%    figure ; plot(xGL,yGL,'ob')
    
    % take out points defining GL that are within the distance CtrlVar.MeshSizeMin of the mesh boundary
    % (improve on this later)
    yMBmin=min(MeshBoundaryCoordinates(:,2)); yMBmax=max(MeshBoundaryCoordinates(:,2));
    ind=(yGL < (yMBmax-CtrlVar.MeshSizeMin)) & (yGL > (yMBmin+CtrlVar.MeshSizeMin));
    
    xGLtemp=xGL; yGLtemp=yGL;
    yGL=yGL(ind) ; xGL=xGL(ind);
    
    
  %  hold on
  %  plot(xGL,yGL,'xr')
    
    
    
    n=size(MeshBoundaryCoordinates,1);
    N=size(xGL,1);
    MeshBoundaryCooWithGLcoo=[MeshBoundaryCoordinates; [xGL(:) yGL(:)]];
    
    
    % ea=[4+(1:(N-1)), 4+N:-1:6 ];
    
    ea=[n+(1:(N-1)), n+N:-1:(n+2) ];
    eb=CircShiftVector(ea,1);
    
    ta=1:n ; tb=CircShiftVector(ta,1);
    edge=[ta' tb'; ea' eb'];
    
    %edge=[1 2 ; 2 3 ; 3 4 ; 4 1 ; [ea' eb'] ];
    face=[];
    
    xGL=xGLtemp; yGL=yGLtemp;  % revert to  unmodified GL points
    
    
    

    return
    
end
