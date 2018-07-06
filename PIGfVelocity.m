function [u,v] = PIGfVelocity(X,Y,Vx,Vy,coordinates,NodeList)
    
    
    Vx(isnan(Vx))=0;
    Vy(isnan(Vy))=0;
    
    u=interp2(X,Y,Vx,coordinates(NodeList,1),coordinates(NodeList,2),'linear');
    v=interp2(X,Y,Vy,coordinates(NodeList,1),coordinates(NodeList,2),'linear');
    
    
    
end
