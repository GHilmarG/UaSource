function MeshBoundaryCoordinates=CreateMeshBoundaryCoordinatesForAntartica(CtrlVar)
    
    load BedmapZeroIceThicknessContour.mat
    
    CtrlVar.GLtension=1e-12; % tension of spline, 1: no smoothing; 0: straight line
    CtrlVar.GLds=CtrlVar.MeshSize;
    % CtrlVar.PlotMesh=1;
    % CtrlVar.fidlog=1;
    % % % take out AP
    % % I=xh0<-1.78e6 & yh0> 5.5e5 ;
    % % xh0=xh0(~I) ; yh0=yh0(~I);
    
    [x,y,nx,ny] = Smooth2dPos(xh0,yh0,CtrlVar);
    
%     % x=flipud(x) ; y=flipud(y) ; nx=flipud(nx) ; ny=flipud(ny) ;
%     figure
%     plot(xh0,yh0,'r.-') ; hold on
%     plot(x,y,'b.-') ;
%     axis equal
%     legend('original','smooth & equidistant')
%     
    
    
    MeshBoundaryCoordinates=[x(:) y(:)];
    %MeshBoundaryCoordinates=flipud(MeshBoundaryCoordinates);

    
end
