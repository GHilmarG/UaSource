function [coordinates10,connectivity10]=tri3to10(coordinates,connectivity)
    % creates 6 node triangles from 3 node tringles
    % reasonably fast
    
    Nnodes=length(coordinates);
    [Nele,nod]=size(connectivity);
    
    %TR = TriRep(connectivity, coordinates(:,1),coordinates(:,2));
    
    % keep it simple, just add an addional midnodes to all edges of all elements, and then get rid of
    % duplicates afterwards
    
    
    nod10=10 ; % number of nodes per element
    
    connectivity10=zeros(Nele,nod10);
    coordinates10=zeros(Nnodes+7*Nele,2);
    
    
    % use same node numbers for corner nodes
    connectivity10(:,1)=connectivity(:,1);
    connectivity10(:,4)=connectivity(:,2);
    connectivity10(:,7)=connectivity(:,3);
	
    coordinates10(connectivity10(:,1),:)=coordinates(connectivity(:,1),:);
    coordinates10(connectivity10(:,4),:)=coordinates(connectivity(:,2),:);
    coordinates10(connectivity10(:,7),:)=coordinates(connectivity(:,3),:);
    
    % create new node labels for the mid nodes
	
	I=Nnodes;
	for Inod=[2 3 5 6 8 9 10];
		newnode=[I+1:I+Nele]';
		I=I+Nele;
		connectivity10(:,Inod)=newnode;
	end
	
	
    % calculate coordinates of mid nodes
    
    coordinates10(connectivity10(:,2),1)=(2*coordinates(connectivity(:,1),1)+coordinates(connectivity(:,2),1))/3;
    coordinates10(connectivity10(:,2),2)=(2*coordinates(connectivity(:,1),2)+coordinates(connectivity(:,2),2))/3;
	
	coordinates10(connectivity10(:,3),1)=(coordinates(connectivity(:,1),1)+2*coordinates(connectivity(:,2),1))/3;
    coordinates10(connectivity10(:,3),2)=(coordinates(connectivity(:,1),2)+2*coordinates(connectivity(:,2),2))/3;
	
	coordinates10(connectivity10(:,5),1)=(2*coordinates(connectivity(:,2),1)+coordinates(connectivity(:,3),1))/3;
    coordinates10(connectivity10(:,5),2)=(2*coordinates(connectivity(:,2),2)+coordinates(connectivity(:,3),2))/3;
	
    coordinates10(connectivity10(:,6),1)=(coordinates(connectivity(:,2),1)+2*coordinates(connectivity(:,3),1))/3;
    coordinates10(connectivity10(:,6),2)=(coordinates(connectivity(:,2),2)+2*coordinates(connectivity(:,3),2))/3;
	
	coordinates10(connectivity10(:,8),1)=(2*coordinates(connectivity(:,3),1)+coordinates(connectivity(:,1),1))/3;
    coordinates10(connectivity10(:,8),2)=(2*coordinates(connectivity(:,3),2)+coordinates(connectivity(:,1),2))/3;
	
    coordinates10(connectivity10(:,9),1)=(coordinates(connectivity(:,3),1)+2*coordinates(connectivity(:,1),1))/3;
    coordinates10(connectivity10(:,9),2)=(coordinates(connectivity(:,3),2)+2*coordinates(connectivity(:,1),2))/3;
	
	coordinates10(connectivity10(:,10),1)=(coordinates(connectivity(:,1),1)+coordinates(connectivity(:,2),1)+coordinates(connectivity(:,3),1))/3;
    coordinates10(connectivity10(:,10),2)=(coordinates(connectivity(:,1),2)+coordinates(connectivity(:,2),2)+coordinates(connectivity(:,3),2))/3;
	
    
    tolerance=1000*eps;
		
    [coordinates10,connectivity10]=RemoveDuplicateNodes(coordinates10,connectivity10,tolerance);
    
    
end

