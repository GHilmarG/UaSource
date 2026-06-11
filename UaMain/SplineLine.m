function [x,y,nx,ny]=SplineLine(x,y,CtrlVar)
    
    
    % Spline approx to (x,y), sampled at equally spaced locations
    % with normals
    
       
    [x,y,I] = Arrange2dPos(x,y);
    [x,y,nx,ny] = Smooth2dPos(x,y,CtrlVar);
    
    if nargout==1
        x.x=x;
        x.y=y;
        x.nx=nx;
        x.ny=ny;
    end
    
end