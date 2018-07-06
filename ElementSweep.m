 function  [coordinates,connectivity] =  ElementSweep(coordinates,connectivity,alpha)
    
   
    [Nele,nod]=size(connectivity);
    xcentre=mean(reshape(coordinates(connectivity,1),Nele,nod),2);
    ycentre=mean(reshape(coordinates(connectivity,2),Nele,nod),2);
    temp=xcentre*cos(alpha)+ycentre*sin(alpha);
    [t,p]=sort(temp);
	
    connectivity=connectivity(p,:);
    
		
	
end

