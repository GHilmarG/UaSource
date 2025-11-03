function ind=CornerNodes(connectivity)
    
    
    [Nele,nod]=size(connectivity);
    
    switch nod
        case 3
            ind=1:Nele;
        case 6
             ind=unique([connectivity(:,1);connectivity(:,3);connectivity(:,5)]);
        case 10
            ind=unique([connectivity(:,1);connectivity(:,4);connectivity(:,7)]);
    end
    
    
end
