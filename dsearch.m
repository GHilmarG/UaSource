function pid=dsearch(x,y,TRI,xi,yi,S,dt)
    
    % simple replacement for dsearch
    % TRI and S are kept as input parameters for compatibility
    % dt is calculated if not provided
    
    fprintf('dsearch replacement \n')
    
    
    if nargin < 7
        dt = DelaunayTri(x(:),y(:));
    end
    
    qrypts=[xi(:) yi(:)];
    pid = nearestNeighbor(dt, qrypts);
    
    pid=pid';
    
end