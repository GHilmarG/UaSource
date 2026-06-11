function [GLdescriptors,dtGL]=GLmorphingInit(CtrlVar,Experiment,coordinates,connectivity,GF,MeshBoundaryCoordinates)
    
    
    %
    % initial step needed for GL morphing.
    % Creates background triangulation based on GL position
    
    % find GL and create GL edges
    [MeshBoundaryCooWithGLcoo,edge]=glLineEdgesFaces(GF,coordinates,connectivity,MeshBoundaryCoordinates,CtrlVar);
    %[coordinates,connectivity]=genmesh2d(Experiment,MeshBoundaryCooWithGLcoo,CtrlVar,edge,face);
    
    
    % create background GL triangulation
    GLConstraints=edge;
    dtGL = DelaunayTri(MeshBoundaryCooWithGLcoo(:,1),MeshBoundaryCooWithGLcoo(:,2),GLConstraints);
    
    % plot background triangulation
    figure ; clf; triplot(dtGL); axis equal;
    %%
   
    GLdescriptors.tri = dtGL.pointLocation(coordinates(:,1), coordinates(:,2));
    GLdescriptors.baryCoords = dtGL.cartToBary(GLdescriptors.tri, coordinates);
    
 
end