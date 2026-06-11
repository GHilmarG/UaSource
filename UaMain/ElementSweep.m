 function  [coordinates,connectivity] =  ElementSweep(coordinates,connectivity,alpha)
    
 % coordinates are used but not changed, not sure why it is here as an output...
   
    [Nele,nod]=size(connectivity);
    xcentre=mean(reshape(coordinates(connectivity,1),Nele,nod),2);
    ycentre=mean(reshape(coordinates(connectivity,2),Nele,nod),2);
    temp=xcentre*cos(alpha)+ycentre*sin(alpha);
    [~,p]=sort(temp);
	
    connectivity=connectivity(p,:);
    
		
	
end

