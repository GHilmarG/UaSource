


function EdgeWork=EdgeWorkIntegral(CtrlVar,MUA,Displacement,Txy,options)

%%
% Calculates the mesh boundary integral:
%
% $$ \int_{\Gamma} \mathbf{s}(\mathbf{n}) \cdot \mathbf{d} \; d \Gamma $$
%
% Although, here it is assumed that 
% 
% $$\mathbf{s}(\mathbf{n})$$
% 
% can be written as
%
% $$T_{xy} \, \mathbf{n}$$
%
% where $T_{xy}$ is a scalar. So what is calculated is:
%
% $$ \int_{\Gamma} ( T_{xy} \, \mathbf{n} ) \cdot \mathbf{d} \; d \Gamma $$
%
% Displacement : nodal variable of displacements across the whole mesh
% 
% Txy          : This is a scalar variable (it is here assumed that the traction is directed normal to the element
%                edges
%
% The approach used here is elementary and is not based on a FE assembly
%
%
% This works for 3, 6 and 10 node elements.
%
% *Output:* 
%
% EdgeWork.Value : An n x m array, where n is the number of edges, and m number of sub-edges. Every triangular element has
% three edges, and 3, 6 and 10 node elements have 1 , 2, and 3 sub-edges, respectively.
%
% EdgeWork.x and EdgeWork.y : x and y coordinates of the (sub) edges.
%
% Example:
%
% EdgeWork=EdgeWorkIntegral(CtrlVar,MUA,Displacement,Txy,Plots=true,Test=true) ;

% 
%%

arguments
    CtrlVar struct
    MUA     struct
    Displacement {mustBeNumeric} 
    Txy {mustBeNumeric} 
    options.Plots logical = false 
    options.Test = false
end


%%

if options.Test

    UserVar=[];
    CtrlVar=Ua2D_DefaultParameters(); %
    CtrlVar.PlotXYscale=1;
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=1;

    % three lines are required in the Ua2D_InitialUserInput.m
    CtrlVar.MeshSizeMax=10;
    CtrlVar.MeshSizeMin=0.1;
    CtrlVar.MeshSize=1;
    CtrlVar.TriNodes=6;

    MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 0];
    MeshBoundaryCoordinates=[-1 -1 ; 1 -1 ; 1 1 ; -1 1 ];
    MeshBoundaryCoordinates=[-2 -1 ; 1 -1 ; 0 1 ];

    CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
    [UserVar,MUA]=genmesh2d(UserVar,CtrlVar);
    Displacement=zeros(MUA.Nnodes,2); % vector quantity
    Displacement(:,1)=1;
    Displacement(:,2)=0;
    Txy=zeros(MUA.Nnodes,1); % scalar quantity
    Txy(:,1)=1;

end

%%

Edges=MUA.Boundary.Edges ; 
[nEdges,nEdgeNodes]=size(Edges); 
nSubEdges=nEdgeNodes-1 ;
Ax=zeros(nEdges,nSubEdges) ;  Ay=zeros(nEdges,nSubEdges) ;
Bx=zeros(nEdges,nSubEdges) ;  By=zeros(nEdges,nSubEdges) ;


for I=1:nSubEdges
    Ax(:,I)=MUA.coordinates(Edges(:,I),1)   ;  Ay(:,I)=MUA.coordinates(Edges(:,I),2);
    Bx(:,I)=MUA.coordinates(Edges(:,I+1),1) ;  By(:,I)=MUA.coordinates(Edges(:,I+1),2);
end

dx=Bx-Ax ; dy=By-Ay;  % edge distances in x and y direction
nx=-dy ; ny=dx ;
lEdge=sqrt(nx.*nx+ny.*ny); % edge length
nx=nx./lEdge ; ny=ny./lEdge;


xEdgeNodalDisplacements=Edges*0;
yEdgeNodalDisplacements=Edges*0;

for I=1:nEdgeNodes

    xEdgeNodalDisplacements(:,I)=Displacement(Edges(:,I),1) ;
    yEdgeNodalDisplacements(:,I)=Displacement(Edges(:,I),2) ;

end

xEdgeTractionNodal=Edges*0;
yEdgeTractionNodal=Edges*0;
% Here I make, the well justified, simplification that the normals to each sub-edge are equal. This is always going to be the
% case because all edges of triangular elements are always straight. 
for I=1:nEdgeNodes
    xEdgeTractionNodal(:,I)=Txy(Edges(:,I)).*nx(:,1) ;
    yEdgeTractionNodal(:,I)=Txy(Edges(:,I)).*ny(:,1) ;
end

EdgesWork=zeros(nEdges,nSubEdges); 


for J=1:nSubEdges

    EdgesWork(:,J)= ...
        (xEdgeTractionNodal(:,J).*xEdgeNodalDisplacements(:,J)  +  xEdgeTractionNodal(:,J+1).*xEdgeNodalDisplacements(:,J+1) + ...
         yEdgeTractionNodal(:,J).*yEdgeNodalDisplacements(:,J)  +  yEdgeTractionNodal(:,J+1).*yEdgeNodalDisplacements(:,J+1) ) ...
        .*lEdge(:,J)/2 ; 

end

xEdge=(Ax+Bx)/2 ;
yEdge=(Ay+By)/2;

EdgeWork.Value=EdgesWork ;
EdgeWork.x=xEdge;
EdgeWork.y=yEdge;


if options.Plots

    %%
    fEW=FindOrCreateFigure("Mesh") ; clf(fEW)
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=false;
    CtrlVar.PlotNodes=true;

    PlotMuaMesh(CtrlVar,MUA)
    hold on

    isNeg=EdgeWork.Value<0;
    scale=100;
    scale=1e-3*max(EdgeWork.x-EdgeWork.y)/max(EdgeWork.Value);
    scatter(EdgeWork.x(isNeg)/CtrlVar.PlotXYscale,EdgeWork.y(isNeg)/CtrlVar.PlotXYscale,-scale*EdgeWork.Value(isNeg),"filled","r")
    scatter(EdgeWork.x(~isNeg)/CtrlVar.PlotXYscale,EdgeWork.y(~isNeg)/CtrlVar.PlotXYscale,scale*EdgeWork.Value(~isNeg)+eps,"filled","b")
    UaPlots(CtrlVar,MUA,[],Displacement,CreateNewFigure=false)
    title("Displacements and boundary work")
    subtitle("work loss (blue)/work gain (red)")
    fprintf(" sum of work %g (should be zero for this test case)\n",sum(EdgeWork.Value(:)))
    %%
end



end


