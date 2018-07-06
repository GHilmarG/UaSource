function [coo6,connectivity6]=tri3to6SlowButWorks(coo,connectivity);
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    [Nnodes,temp]=size(coo);
    [Nele,nod]=size(connectivity);
    
    % keep it simple, just add an addional midnodes to all edges of all elementa, and then get rid of
    % duplicates afterwards
    
    
    nod6=6 ; % number of nodes per element
    
    connectivity6=zeros(Nele,nod6);
    coo6=zeros(Nnodes+3*Nele,2);
    
    
    % use same node numbers for corner nodes
    connectivity6(:,1)=connectivity(:,1);
    connectivity6(:,3)=connectivity(:,2);
    connectivity6(:,5)=connectivity(:,3);
    coo6(connectivity6(:,1),:)=coo(connectivity(:,1),:);
    coo6(connectivity6(:,3),:)=coo(connectivity(:,2),:);
    coo6(connectivity6(:,5),:)=coo(connectivity(:,3),:);
    
    % create new node labels for the mid nodes
    
    newnode=Nnodes;

    for Iele=1:Nele
        for Inod=2:2:6
            newnode=newnode+1;
            connectivity6(Iele,Inod)=newnode;
        end
    end
    
    % calculate coordinates of mid nodes
    
    coo6(connectivity6(:,2),1)=(coo(connectivity(:,1),1)+coo(connectivity(:,2),1))/2;
    coo6(connectivity6(:,4),1)=(coo(connectivity(:,2),1)+coo(connectivity(:,3),1))/2;
    coo6(connectivity6(:,6),1)=(coo(connectivity(:,3),1)+coo(connectivity(:,1),1))/2;
    coo6(connectivity6(:,2),2)=(coo(connectivity(:,1),2)+coo(connectivity(:,2),2))/2;
    coo6(connectivity6(:,4),2)=(coo(connectivity(:,2),2)+coo(connectivity(:,3),2))/2;
    coo6(connectivity6(:,6),2)=(coo(connectivity(:,3),2)+coo(connectivity(:,1),2))/2;
    
    % now go through all the new mid nodes, see if a node is a duplicate and if so eliminate
    
    lastnode=Nnodes+3; % don't check first element
    
    for Iele=2:Nele
        
        for midnode=2:2:6;
%             disp(' ');
%             disp('--------------------------------------------------------------')
%             disp([' lastnode : ',num2str(lastnode)])
            dist=zeros(3,1) ; ind=zeros(3,1);
            for J=1:3
                [dist(J),ind(J)]=min(abs(coo6(connectivity6(Iele,midnode),1)-coo6(connectivity6(1:Iele-1,J*2),1))+...
                    abs(coo6(connectivity6(Iele,midnode),2)-coo6(connectivity6(1:Iele-1,J*2),2)));
            end
            
            [Dist,Ind]=min(dist) ;
            
            %disp([' dist : ',num2str(dist(1)),',',num2str(dist(2)),',',num2str(dist(3))])
            %disp([' Dist : ',num2str(Dist)])
            %disp([' Ind : ',num2str(Ind)])
            
            if  Dist < 1000 *eps
                
                node=connectivity6(ind(Ind),Ind*2);
                %disp([' duplicate, replaced by : ',num2str(node)])
                
                
                
            else
                
                node=lastnode+1; lastnode=node;
                %disp([' new node ',num2str(node)])
            end
            
            %connectivity6(Iele,:)
            
            coo6(node,:)=coo6(connectivity6(Iele,midnode),:);
            connectivity6(Iele,midnode)=node;

            %connectivity6(Iele,:)
            %disp([' coordinates : ',num2str(coo6(node,1)),',',num2str(coo6(node,2))])
            
        end
    end
    
    

coo6=coo6(1:lastnode,:);    

end

