function connectivityOut=FlipElements(connectivityIn,I)
    
    %
    % connectivityOut=FlipElements(connectivityIn,I)
    % Reverses the connectvity of elements, ie flips them inside-out.
    % If the optional list I is given then only those elements will be flipped.
    % The list `I' can be either a logical list or a list of elements.
    
    [Nele,nod]=size(connectivityIn);
    connectivityOut=connectivityIn;
    
    if nargin==1 
        I=true(size(connectivityIn,1),1);
    end
    
    switch nod
        
        case 3
            
            connectivityOut(I,1)=connectivityIn(I,1);
            connectivityOut(I,2)=connectivityIn(I,3);
            connectivityOut(I,3)=connectivityIn(I,2);
            
        case 6
            
            connectivityOut(I,1)=connectivityIn(I,1);
            connectivityOut(I,2)=connectivityIn(I,6);
            connectivityOut(I,3)=connectivityIn(I,5);
            connectivityOut(I,4)=connectivityIn(I,4);
            connectivityOut(I,5)=connectivityIn(I,3);
            connectivityOut(I,6)=connectivityIn(I,2);
            
            
        case 10
            
            connectivityOut(I,1)=connectivityIn(I,1);
            connectivityOut(I,2)=connectivityIn(I,9);
            connectivityOut(I,3)=connectivityIn(I,8);
            connectivityOut(I,4)=connectivityIn(I,7);
            connectivityOut(I,5)=connectivityIn(I,6);
            connectivityOut(I,6)=connectivityIn(I,5);
            connectivityOut(I,7)=connectivityIn(I,4);
            connectivityOut(I,8)=connectivityIn(I,3);
            connectivityOut(I,9)=connectivityIn(I,2);
            connectivityOut(I,10)=connectivityIn(I,10);
            
    end
    
    
end


