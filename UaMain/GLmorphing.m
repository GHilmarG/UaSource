function [coordinates,xGLmorphing,yGLmorphing]=GLmorphing(UserVar,CtrlVar,coordinates,connectivity,GF,MeshBoundaryCoordinates,GLdescriptors,dtGL)
    
    
    
    
    [~,~,xGL,yGL]=GLgeometry(connectivity,coordinates,GF,CtrlVar);
   
    % PlotFEmesh(coordinates,connectivity,CtrlVar);
    % hold on ; plot(xGL,yGL,'r','LineWidth',2);
    
    % If a GL is close but not at the margins of the domain shift it onto the margin to ensure
    % that the nodes along the margins move with the GL 
    %
      
        
    yMBmin=min(MeshBoundaryCoordinates(:,2)); yMBmax=max(MeshBoundaryCoordinates(:,2));
    [yGLmin,imin]=min(yGL);
    [yGLmax,imax]=max(yGL);
  
    if yGLmin< (yMBmin+CtrlVar.MeshSizeMin) ; yGL(imin)=yMBmin; end
    if yGLmax> (yMBmax-CtrlVar.MeshSizeMin) ; yGL(imax)=yMBmax; end
       
    xGLmorphing=xGL ; yGLmorphing=yGL;
    
    MeshBoundaryCooWithGLcoo=[MeshBoundaryCoordinates; [xGL(:) yGL(:)]];
    

    trGL = TriRep(dtGL(:,:),MeshBoundaryCooWithGLcoo(:,1),MeshBoundaryCooWithGLcoo(:,2));
  %  triplot(trGL)
    

    coordinates = trGL.baryToCart(GLdescriptors.tri, GLdescriptors.baryCoords);

   % CtrlVar.MeshColor='r';
   % PlotFEmesh(coordinates,connectivity,CtrlVar);
   % hold on ; plot(xGL,yGL,'g','LineWidth',2);
    
%% temporary fix
   ymax=max(MeshBoundaryCoordinates(:,2)); ymin=min(MeshBoundaryCoordinates(:,2));
   
   I=abs(coordinates(:,2)-ymax)<1e-8 ; coordinates(I,2)=ymax; 
   I=abs(coordinates(:,2)-ymin)<1e-8 ; coordinates(I,2)=ymin; 
   
   
   
    
end