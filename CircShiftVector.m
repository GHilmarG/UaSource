function vector = CircShiftVector(vector,ishift)
  
    % circular shift of a vector
            
     ind=mod(ishift+(0:(numel(vector)-1)),numel(vector))+1;
     vector=vector(ind);
         
    
end

