function [BoundaryNodes,BoundaryEdgeCornerNodes,BoundaryFreeElements,BoundaryEdges,BoundaryEdge]=FindBoundaryNodes(connectivity,coordinates)
    
    %  old version, now use FindBoundary.m instead
    % [BoundaryNodes,BoundaryEdgeCornerNodes,BoundaryFreeElements]=FindBoundaryNodes(connectivity,coordinates)
    %
    %           BoundaryNodes : list of all nodes on the boundary, ie not only corner nodes (it is not a linked list).
    % BoundaryEdgeCornerNodes : a linked array of the corner nodes of free edges.
    
    %    BoundaryFreeElements : a list of elements with free edges, each edge listed once, ie if an
    %                           element has more than one free edge, it is listed more than once.
    %           BoundaryEdges : all free boundary edges of the FE mesh, each edged defined by the two corner nodes
    %            BoundaryEdge : cell array with three elements listing the relative node numbers for each edge depending on the type of element.
    %
    % The boundary can be plotted using PlotBoundary(Boundary,connectivity,coordinates,CtrlVar)
    %
        
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    
    if Nele==0;
        BoundaryNodes=[];
        BoundaryEdgeCornerNodes=[];
        BoundaryFreeElements=[];
        BoundaryEdges=[];
        BoundaryEdge=[];
        return
    end
    
    switch nod
        case 3
            con3=connectivity;
            BoundaryEdge{1}=[1 2] ;  BoundaryEdge{2}=[2 3] ;  BoundaryEdge{3}=[3 1] ;
        case 6
            con3=connectivity(:,[1 3 5]);
            BoundaryEdge{1}=[1 2 3] ;  BoundaryEdge{2}=[3 4 5] ;  BoundaryEdge{3}=[5 6 1] ;
        case 10
            con3=connectivity(:,[1 4 7]);
            BoundaryEdge{1}=[1 2 3 4] ;  BoundaryEdge{2}=[4 5 6 7] ;  BoundaryEdge{3}=[7 8 9 1] ;
        otherwise
            error(' case not recognized')
    end
    
    warning('off','MATLAB:TriRep:PtsNotInTriWarnId')
    TR = TriRep(con3,coordinates(:,1),coordinates(:,2));
    
    BoundaryEdges=freeBoundary(TR) ; % misses the interor nodes for higher order tri
    %BoundaryEdgeCornerNodes=unique(BoundaryEdges(:));
    BoundaryEdgeCornerNodes=BoundaryEdges(:,1);
    BoundaryFreeElements=cell2mat(edgeAttachments(TR,BoundaryEdges));  % OK of every type
    % edgeAttachments returns a cell array, but because each boundary edge only belongs to one element
    % each element of that cell array has the same number of elements, or just 1
    nBoundaryEdges=size(BoundaryEdges,1);
    nBoundaryElements=size(BoundaryEdges,1);
    
    edge1=logical(prod(double(con3(BoundaryFreeElements,[1 2]) == BoundaryEdges),2));
    edge2=logical(prod(double(con3(BoundaryFreeElements,[2 3]) == BoundaryEdges),2));
    edge3=logical(prod(double(con3(BoundaryFreeElements,[3 1]) == BoundaryEdges),2));
    
    % I now have to add the boundary nodes for tri6 and tri10
    
    if nod==6
        temp=zeros(nBoundaryEdges,3);
        temp(edge1,:)=connectivity(BoundaryFreeElements(edge1),[1 2 3]);
        temp(edge2,:)=connectivity(BoundaryFreeElements(edge2),[3 4 5]);
        temp(edge3,:)=connectivity(BoundaryFreeElements(edge3),[5 6 1]);
        BoundaryEdges=temp;
    elseif nod==10
        temp=zeros(nBoundaryEdges,4);
        temp(edge1,:)=connectivity(BoundaryFreeElements(edge1),[1 2 3 4]);
        temp(edge2,:)=connectivity(BoundaryFreeElements(edge2),[4 5 6 7]);
        temp(edge3,:)=connectivity(BoundaryFreeElements(edge3),[7 8 9 1]);
        BoundaryEdges=temp;
    end
    
    
    BoundaryNodes=unique(BoundaryEdges);
    
    %Boundary.Edges1Ptrs=edge1;
    %Boundary.Edges2Ptrs=edge2;
    %Boundary.Edges3Ptrs=edge3;
    %Boundary.Edges=BoundaryEdges;
    %Boundary.Elements=BoundaryElements;
    %InteriorNodes=setdiff(1:Nnodes,BoundaryNodes);
    
    
end


