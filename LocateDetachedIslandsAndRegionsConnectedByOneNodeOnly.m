




function [Islands]=LocateDetachedIslandsAndRegionsConnectedByOneNodeOnly(CtrlVar,MUA)





% This works for all element types, but the `implementation' for 6 and 10 nodes is simply based on creating a new 3-node mesh
% with the same element numbering
if CtrlVar.TriNodes~=3
    CtrlVar.TriNodes=3 ;
    CtrlVar.CalcMUA_Derivatives=false;
    CtrlVar.MUA.MassMatrix=false;
    CtrlVar.MUA.DecomposeMassMatrix=false;
    CtrlVar.MUA.StiffnessMatrix=false;
    CtrlVar.FindMUA_Boundary=true;
    MUA=UpdateMUA(CtrlVar,MUA) ;
end


TRI=MUA.connectivity;
G=graph(TRI,TRI(:,[2 3 1]));



FindOrCreateFigure("Graph") ; plot(G)


% calculate the connected components of the graph
[bin,binsize]=conncomp(G) ;
% size of bin is equal length to the number of nodes
% bin(i) is the bin number to which node i belongs

idx = binsize(bin) == max(binsize);
nodes=find(~idx);   % all nodes not belonging to the subgraph with the largest number of nodes

ElementsToBeDeleted=AllElementsContainingGivenNodes(MUA.connectivity,nodes) ;

Islands.Free=ElementsToBeDeleted ;




% 2) Get rid of elements only connected by one node to the main part of the mesh


% Find boundary node that appears more than once along the boundary
V=MUA.Boundary.EdgeCornerNodes ; sV=sort(V) ; BoundaryNodeDuplicates=unique(sV(diff(sV)==0)) ;

% Replace each occurence of those duplicates with new node numbers


TempConnectivity=MUA.connectivity ;
NodeMax=MUA.Nnodes;
for I=1:numel(BoundaryNodeDuplicates)

    ind=find(TempConnectivity==BoundaryNodeDuplicates(I)) ;
    ind=ind(2:end) ;  % do not relabel the first instance of the node
    n=numel(find(ind));

    NewNodeLables=(1:n)+NodeMax ;
    TempConnectivity(ind)=NewNodeLables;
    NodeMax=NodeMax+n;

    % I only need this for plotting purposes
    %for J=1:numel(NewNodeLables)
    %    MUA.coordinates(NewNodeLables(J),:)=MUA.coordinates(BoundaryNodeDuplicates(I),:);
    %end

end

% Now the connectivity have changed, but the number of elements is the same



if CtrlVar.TriNodes==6
    TRI=TempConnectivity(:,[1 3 5]);
elseif CtrlVar.TriNodes==10
    TRI=TempConnectivity(:,[1 4 7]);
else
    TRI=TempConnectivity;
end

G=graph(TRI,TRI(:,[2 3 1]));


% FindOrCreateFigure("Graph 2") ; plot(G)

[bin,binsize]=conncomp(G) ;


idx = binsize(bin) == max(binsize);
nodes=find(~idx);   % all nodes not belonging to the subgraph with the largest number of nodes


FurtherElementsToBeDeleted=AllElementsContainingGivenNodes(TempConnectivity,nodes) ;

Islands.OneNode=FurtherElementsToBeDeleted ;

% Since the element labels and number of elements is the same, I can now simply use the original MUA

%%



end

%%