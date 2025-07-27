
%%
UserVar=[];
CtrlVar=Ua2D_DefaultParameters(); %
CtrlVar.PlotXYscale=1; 
 CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=1;
% Note; When creating this mesh using Ua, only the following 
% three lines are required in the Ua2D_InitialUserInput.m
CtrlVar.MeshSizeMax=10; 
CtrlVar.MeshSizeMin=0.1;
CtrlVar.MeshSize=1;
CtrlVar.TriNodes=10; 

MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 0];
MeshBoundaryCoordinates=[-1 -1 ; 1 -1 ; 1 1 ; -1 1 ];
MeshBoundaryCoordinates=[-2 -1 ; 1 -1 ; 0 1 ];

CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
% Now generate mesh (When using Ua this is done internally, no such call
% then needed).


[UserVar,MUA]=genmesh2d(UserVar,CtrlVar);



Displacement=zeros(MUA.Nnodes,2); % vector quantity 
Displacement(:,1)=0;
Displacement(:,2)=1;
Txy=zeros(MUA.Nnodes,1); % scalar quantity 
Txy(:,1)=1;


 EdgeWork=EdgeWorkIntegral(CtrlVar,MUA,Displacement,Txy,Plots=true) ;

 


 %%
