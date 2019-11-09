function h = hfun(x,y,MeshSize,MeshSizeMin,MeshSizeMax)
    
    
    if numel(MeshSize) >1
        
        
        % F = TriScatteredInterp(MeshSize(:,1),MeshSize(:,2),MeshSize(:,3),'natural');
        F = scatteredInterpolant(MeshSize(:,1),MeshSize(:,2),MeshSize(:,3),'natural');
        F.Method = 'natural';
        h = F(x,y);
        
    else
        h =  MeshSize+abs(x)*0.;
    end
    
    if nargin > 3
        h(h<MeshSizeMin)=MeshSizeMin ; h(h>MeshSizeMax)=MeshSizeMax;
    end
    
    
end

