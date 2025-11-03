function  [dtGL,GLdescriptors,xGLbackground,yGLbackground]=...
        CreateBackgroundGLmesh(coordinates,MeshBoundaryCoordinates,xGLmesh,yGLmesh,CtrlVar)
    
    %   Detailed explanation goes here
    
    % if the GL point closest to the boundary of the mesh, is within the distance CtrlVar.MeshSizeMin from the boundary
    % move it so that it is exactly on the boundary.  This is done to ensure that the background triangulation
    % in that case has a point on the boundary that moves with the GL
    % (currently just considering at simple box shaped mesh domain, improve on this later)
    
    yMBmin=min(MeshBoundaryCoordinates(:,2)); yMBmax=max(MeshBoundaryCoordinates(:,2));
    [yGLmin,imin]=min(yGLmesh);
    [yGLmax,imax]=max(yGLmesh);
  
    % if GL point close to boundary, put exactly on boundary
    if yGLmin< (yMBmin+CtrlVar.MeshSizeMin) ; yGLmesh(imin)=yMBmin; end
    if yGLmax> (yMBmax-CtrlVar.MeshSizeMin) ; yGLmesh(imax)=yMBmax; end
    
    xGLbackground=xGLmesh ; yGLbackground=yGLmesh;
    
    n=size(MeshBoundaryCoordinates,1);
    
    MeshBoundaryCooWithGLcoo=[MeshBoundaryCoordinates; [xGLmesh(:) yGLmesh(:)]];
    
        
    
    ta=1:n ; tb=CircShiftVector(ta,1);
        
    
    
    GLConstraints=[ta' tb'];
    dtGL = DelaunayTri(MeshBoundaryCooWithGLcoo(:,1),MeshBoundaryCooWithGLcoo(:,2),GLConstraints);
    %dtGL = DelaunayTri(MeshBoundaryCooWithGLcoo(:,1),MeshBoundaryCooWithGLcoo(:,2));
    
    
    % plot background triangulation
    %figure ; clf; triplot(dtGL); axis equal;
    
    
    GLdescriptors.tri = dtGL.pointLocation(coordinates(:,1), coordinates(:,2));
    GLdescriptors.baryCoords = dtGL.cartToBary(GLdescriptors.tri, coordinates);
    
    
    
    
end

