function A=TriArea(x,y)
    
 % calculates area of multiple number of triangles   
 % 
 %  X=[x1 x2 x3 ]; Y=[y1 y2 y3] where each component is a vector
 %
 
 
    A=0.5*abs((x(:,1)-x(:,3)).*(y(:,2)-y(:,1)) -(x(:,1)-x(:,2)).*(y(:,3)-y(:,1)));
    
end
