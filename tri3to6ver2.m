function [coordinates6,connectivity6]=tri3to6ver2(coordinates,connectivity);
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    [Nnodes,temp]=size(coordinates);
    [Nele,nod]=size(connectivity);
    
    TR = TriRep(connectivity, coordinates(:,1),coordinates(:,2));
    
    % keep it simple, just add an addional midnodes to all edges of all elementa, and then get rid of
    % duplicates afterwards
    
    
    nod6=6 ; % number of nodes per element
    
    connectivity6=zeros(Nele,nod6);
    coordinates6=zeros(Nnodes+3*Nele,2);
    
    
    % use same node numbers for corner nodes
    connectivity6(:,1)=connectivity(:,1);
    connectivity6(:,3)=connectivity(:,2);
    connectivity6(:,5)=connectivity(:,3);
    coordinates6(connectivity6(:,1),:)=coordinates(connectivity(:,1),:);
    coordinates6(connectivity6(:,3),:)=coordinates(connectivity(:,2),:);
    coordinates6(connectivity6(:,5),:)=coordinates(connectivity(:,3),:);
    
    % create new node labels for the mid nodes
    
    newnode=Nnodes;

    for Iele=1:Nele
        for Inod=2:2:6
            newnode=newnode+1;
            connectivity6(Iele,Inod)=newnode;
        end
    end
    
    % calculate coordinates of mid nodes
    
    coordinates6(connectivity6(:,2),1)=(coordinates(connectivity(:,1),1)+coordinates(connectivity(:,2),1))/2;
    coordinates6(connectivity6(:,4),1)=(coordinates(connectivity(:,2),1)+coordinates(connectivity(:,3),1))/2;
    coordinates6(connectivity6(:,6),1)=(coordinates(connectivity(:,3),1)+coordinates(connectivity(:,1),1))/2;
    coordinates6(connectivity6(:,2),2)=(coordinates(connectivity(:,1),2)+coordinates(connectivity(:,2),2))/2;
    coordinates6(connectivity6(:,4),2)=(coordinates(connectivity(:,2),2)+coordinates(connectivity(:,3),2))/2;
    coordinates6(connectivity6(:,6),2)=(coordinates(connectivity(:,3),2)+coordinates(connectivity(:,1),2))/2;
    
    
    tolerance=1;
    [coordinates6,connectivity6]=RemoveDuplicateNodes(coordinates6,connectivity6,tolerance);
    
    
end

