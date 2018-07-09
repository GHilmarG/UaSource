function [coordinates,connectivity]=NodalSweep(coordinates,connectivity,alpha)
        
    temp=coordinates(:,1)*cos(alpha)+coordinates(:,2)*sin(alpha);
    [t,p]=sort(temp);
    
    p2=p*0; p2(p)=1:length(p2);
    coordinates(:,2)=coordinates(p,2);
    coordinates(:,1)=coordinates(p,1);
    connectivity=p2(connectivity);
            
    
end

