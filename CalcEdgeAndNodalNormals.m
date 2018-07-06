function [nx,ny,xn,yn,Nx,Ny] = CalcEdgeAndNodalNormals(connectivity,coordinates,Edges)

%%
% [nx,ny,xn,yn] = EdgeNormals(connectivity,coordinates,Edges)
% calculates edge and (optionally) nodal normals
%
%  nx, ny : x, y components of the edge normals
%  xn, yn : x, y coordinates of the edge normals (centre of the edge)
%
%  Nx,Ny  : nodal normals
%
%
%
% Edges are defined by listing the nodes along each edge in the array `Edges'.
% The format is Edges=[node1 node2 ; node1 node2 ; etc ] for a 3-node element
%               Edges=[node1 node2 node3 ; node1 node2 node3 ; etc ] for a 6-node element
%               Edges=[node1 node2 node3 node4 ; node1 node2 node3  node4 ; etc ] for a 10-node element
%
% For 6 and 10 node elements, each edge therefore consists of 2 and 3, respectivily, sub-edges.
%
% MUA.Boundary.Edges is in this format
%
% Note that since for higher order elements MUA.Boundary.Edges returns edges and sub-edges
% nx, ny, xn, and yn are returned as #Edges x #NumSubEdges
%
%
% Nx and Ny are optional outputs. These are arrays with equal elements as the total number of nodes
% in mesh. For nodes not in Edges, Nx and Ny are not calculated and NaN is returned
% Nodal normals are defined as the average of the edge normals on both the sides of the node.
% Although nodal normals can be calculated for all nodes, it generally only makes sense to calculate
% them for (free) boundary nodes.
%
% Examples:
%
% To calculate and plot normals to all free edges:
% [nx,ny,xn,yn] = CalcEdgeAndNodalNormals(MUA.connectivity,MUA.coordinates,MUA.Boundary.Edges);
% figure ; QuiverColorGHG(nx,ny,xn,yn);
%
% To calculate and plot nodal normals for boundary nodes:
% [nx,ny,xn,yn,Nx,Ny] = CalcEdgeAndNodalNormals(MUA.connectivity,MUA.coordinates,MUA.Boundary.Edges);
% QuiverColorGHG(MUA.coordinates(MUA.Boundary.Nodes,1),MUA.coordinates(MUA.Boundary.Nodes,2),...
%    Nx(MUA.Boundary.Nodes),Ny(MUA.Boundary.Nodes),CtrlVar);
%%

[nEdges,nod]=size(Edges); nSubEdges=nod-1 ;
Ax=zeros(nEdges,nSubEdges) ;  Ay=zeros(nEdges,nSubEdges) ;
Bx=zeros(nEdges,nSubEdges) ;  By=zeros(nEdges,nSubEdges) ;



for I=1:nSubEdges
    Ax(:,I)=coordinates(Edges(:,I),1);  Ay(:,I)=coordinates(Edges(:,I),2);
    Bx(:,I)=coordinates(Edges(:,I+1),1);  By(:,I)=coordinates(Edges(:,I+1),2);
end

dx=Bx-Ax ; dy=By-Ay;
nx=-dy ; ny=dx ;
l=sqrt(nx.*nx+ny.*ny);
nx=nx./l ; ny=ny./l ;

xn=(Ax+Bx)/2; yn=(Ay+By)/2;


if nargout> 4  % calculate normals for nodes
    Nx=zeros(length(coordinates),1); Ny=zeros(length(coordinates),1); % all nodes!
    Nx(Edges(:,1))=nx(:,1) ; Ny(Edges(:,1))=ny(:,1) ;
    
    
    Nx(Edges(:,nod))=Nx(Edges(:,nod))+ nx(:,nSubEdges) ; % adding up normals
    Ny(Edges(:,nod))=Ny(Edges(:,nod))+ ny(:,nSubEdges) ;
    
    for I=2:nSubEdges
        
        Nx(Edges(:,I))=Nx(Edges(:,I))+ ( nx(:,I-1) + nx(:,I))/2 ;
        Ny(Edges(:,I))=Ny(Edges(:,I))+ ( ny(:,I-1) + ny(:,I))/2 ;
        
    end
    
    l=sqrt(Nx.*Nx+Ny.*Ny);
    Nx=Nx./l ; Ny=Ny./l ;
    
end

end