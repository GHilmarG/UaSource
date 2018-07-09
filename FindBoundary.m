function [Boundary,TR]=FindBoundary(connectivity,coordinates)

%
% [Boundary.Nodes,Boundary.EdgeCornerNodes,Boundary.Free,Elements]=FindBoundary(connectivity,coordinates)
%
%           Boundary.Nodes : list of all nodes on the boundary, ie not only corner nodes (it is not a linked list).
% Boundary.EdgeCornerNodes : a linked array of the corner nodes of free edges.
%
%    Boundary.FreeElements : a list of elements with free edges, each edge listed once, ie if an
%                            element has more than one free edge, it is listed more than once.
%           Boundary.Edges : all free boundary edges of the FE mesh, each edge defined by the two corner nodes
%            Boundary.Edge : cell array with three elements listing the relative node numbers for each edge depending on the type of element.
%
% The boundary can be plotted using PlotBoundary(Boundary,connectivity,coordinates,CtrlVar)
%
% Note that the first and last of Boundary.EdgeCornerNodes are not the same nodes
% i.e. the loop does note close. If doing inside-out tests, then this loops must be closed first
%


[Nele,nod]=size(connectivity);

if Nele==0
    Boundary.Nodes=[];
    Boundary.EdgeCornerNodes=[];
    Boundary.FreeElements=[];
    Boundary.Edges=[];
    Boundary.Edge=[];
    Boundary.x=[]; Boundary.y=[];
    TR=[];
    return
end


switch nod
    case 3
        con3=connectivity;
        Boundary.Edge{1}=[1 2] ;  Boundary.Edge{2}=[2 3] ;  Boundary.Edge{3}=[3 1] ;
    case 6
        con3=connectivity(:,[1 3 5]);
        Boundary.Edge{1}=[1 2 3] ;  Boundary.Edge{2}=[3 4 5] ;  Boundary.Edge{3}=[5 6 1] ;
    case 10
        con3=connectivity(:,[1 4 7]);
        Boundary.Edge{1}=[1 2 3 4] ;  Boundary.Edge{2}=[4 5 6 7] ;  Boundary.Edge{3}=[7 8 9 1] ;
    otherwise
        error(' case not recognized')
end



%TR=CreateFEmeshCornerPointTriangulation(connectivity,coordinates);
TR=CreateFEmeshTriRep(connectivity,coordinates);

Boundary.Edges=freeBoundary(TR) ; % misses the interor nodes for higher order tri
%Boundary.EdgeCornerNodes=unique(Boundary.Edges(:));
Boundary.EdgeCornerNodes=Boundary.Edges(:,1);
Boundary.FreeElements=cell2mat(edgeAttachments(TR,Boundary.Edges));  % OK of every type
% edgeAttachments returns a cell array, but because each boundary edge only belongs to one element
% each element of that cell array has the same number of elements, or just 1
nBoundary.Edges=size(Boundary.Edges,1);
nBoundaryElements=size(Boundary.Edges,1);

edge1=logical(prod(double(con3(Boundary.FreeElements,[1 2]) == Boundary.Edges),2));
edge2=logical(prod(double(con3(Boundary.FreeElements,[2 3]) == Boundary.Edges),2));
edge3=logical(prod(double(con3(Boundary.FreeElements,[3 1]) == Boundary.Edges),2));

% I now have to add the boundary nodes for tri6 and tri10

if nod==6
    temp=zeros(nBoundary.Edges,3);
    temp(edge1,:)=connectivity(Boundary.FreeElements(edge1),[1 2 3]);
    temp(edge2,:)=connectivity(Boundary.FreeElements(edge2),[3 4 5]);
    temp(edge3,:)=connectivity(Boundary.FreeElements(edge3),[5 6 1]);
    Boundary.Edges=temp;
elseif nod==10
    temp=zeros(nBoundary.Edges,4);
    temp(edge1,:)=connectivity(Boundary.FreeElements(edge1),[1 2 3 4]);
    temp(edge2,:)=connectivity(Boundary.FreeElements(edge2),[4 5 6 7]);
    temp(edge3,:)=connectivity(Boundary.FreeElements(edge3),[7 8 9 1]);
    Boundary.Edges=temp;
    
    
end


Boundary.Nodes=unique(Boundary.Edges);


for II=1:3

    Boundary.Elements{II}=find(logical(prod(double(ismember(connectivity(:,Boundary.Edge{II}),Boundary.Nodes)')))); %
end

xa=coordinates(Boundary.Edges(:,1),1); xb=coordinates(Boundary.Edges(:,end),1);
ya=coordinates(Boundary.Edges(:,1),2); yb=coordinates(Boundary.Edges(:,end),2);

[Boundary.x,Boundary.y]=LineUpEdges2([],xa,xb,ya,yb);



end


