function dtcritical=CalcCFLdt2D(DTxy,coordinates,connectivity,u,v);
    
    
    % just a rough method of estimating maximum time based on CFL
    % must improve on this later
    
    
    [Nele,nod]=size(connectivity);
    
    x=coordinates(:,1); y=coordinates(:,2);
    
    area=(max(x)-min(x))*(max(y)-min(y));
    
    dl=sqrt(area/Nele);
    
    dtcritical=dl/max(sqrt(u.*u+v.*v));
    
end

