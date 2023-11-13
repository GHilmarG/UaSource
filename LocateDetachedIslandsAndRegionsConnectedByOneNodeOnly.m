




function [Islands]=LocateDetachedIslandsAndRegionsConnectedByOneNodeOnly(CtrlVar,MUA)


if ~isfield(CtrlVar,"LocateDetachedIslandsAndRegionsConnectedByOneNodeOnly")

  CtrlVar.LocateDetachedIslandsAndRegionsConnectedByOneNodeOnly="-Islands-OneNodeOrLessConnections-" ;

end

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

if ~contains(CtrlVar.LocateDetachedIslandsAndRegionsConnectedByOneNodeOnly,"-OneNodeOrLessConnections-")

    TRI=MUA.connectivity;
    G=graph(TRI,TRI(:,[2 3 1]));



    % FindOrCreateFigure("Graph") ; plot(G)


    % calculate the connected components of the graph
    [bin,binsize]=conncomp(G) ;
    % size of bin is equal length to the number of nodes
    % bin(i) is the bin number to which node i belongs

    idx = binsize(bin) == max(binsize);
    nodes=find(~idx);   % all nodes not belonging to the subgraph with the largest number of nodes

    ElementsToBeDeleted=AllElementsContainingGivenNodes(MUA.connectivity,nodes) ;

    Islands.Free=ElementsToBeDeleted ;
    Islands.OneNode=[];

else


    % 2) Get rid of elements only connected by one node to the main part of the mesh
    %
    % The key idea is to identify boundary nodes that appear twice or more. At those node the freeBoundary of the triangulation
    % `crosses over'.  Each such node is then split up in several new nodes by giving new nodal labels to the duplicates (keeping the
    % first instance unchanged). This changes the topology of the mesh. Only the connectivity is affected, not the
    % coordinates, which play no role. The effect is that the mesh, as a graph, is now detached at this node. A graph is created
    % for this new connectivity and islands detected in the usual way.
    %
    % Note that this will always include all the original islands as well. So if we want to get rid of all islands and all
    % small regions connected by just one node to the remaining part, we only need to to this step.
    %
    % Within MUA we have the field
    %
    %   MUA.Boundary.EdgeCornerNodes
    %
    % which is calculagted as
    %
    %   TR=CreateFEmeshTriRep(connectivity,coordinates);
    %   Boundary.Edges=freeBoundary(TR) ;
    %   Boundary.EdgeCornerNodes=Boundary.Edges(:,1);
    %
    % Find boundary node that appears more than once along the boundary edges

    % Find duplicates

    if isempty(MUA.Boundary)
        %[MUA.Boundary,MUA.TR]=FindBoundary(MUA.connectivity,MUA.coordinates);
        MUA.TR=CreateFEmeshTriRep(MUA.connectivity,MUA.coordinates);
        MUA.Boundary.Edges=freeBoundary(MUA.TR) ; % misses the interor nodes for higher order tri
        MUA.Boundary.EdgeCornerNodes=MUA.Boundary.Edges(:,1);
    end

    V=MUA.Boundary.EdgeCornerNodes ;
    sV=sort(V) ;
    BoundaryNodeDuplicates=unique(sV(diff(sV)==0)) ;

    % Replace each occurence of those duplicates with new node numbers


    TempConnectivity=MUA.connectivity ;
    NodeMax=MUA.Nnodes;
    for I=1:numel(BoundaryNodeDuplicates)

        ind=find(TempConnectivity==BoundaryNodeDuplicates(I)) ;
        ind=ind(2:end) ;  % do not relabel the first instance of the node
        n=numel(ind);

        NewNodeLables=(1:n)+NodeMax ;
        TempConnectivity(ind)=NewNodeLables;
        NodeMax=NodeMax+n;

        % I only need this for plotting purposes
        %for J=1:numel(NewNodeLables)
        %    MUA.coordinates(NewNodeLables(J),:)=MUA.coordinates(BoundaryNodeDuplicates(I),:);
        %end

    end

    % Now the connectivity have changed, but the number of elements is the same


    G=graph(TempConnectivity,TempConnectivity(:,[2 3 1]));


    % FindOrCreateFigure("Graph 2") ; plot(G)

    [bin,binsize]=conncomp(G) ;


    idx = binsize(bin) == max(binsize);
    nodes=find(~idx);   % all nodes not belonging to the subgraph with the largest number of nodes


    FurtherElementsToBeDeleted=AllElementsContainingGivenNodes(TempConnectivity,nodes) ;

    Islands.Free=[]; 
    Islands.OneNode=FurtherElementsToBeDeleted ;

    % Since the element labels and number of elements is the same, I can now simply use the original MUA

    %%

end

end

%%