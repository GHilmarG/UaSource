%%
CtrlVar=Ua2D_DefaultParameters();

% Grid 1
CtrlVar.MeshSizeMax=0.1;
CtrlVar.MeshSizeMin=0.1;
MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1];


MUA1=genmesh2d(CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA1.coordinates,MUA1.connectivity)

h1=2*abs(MUA1.coordinates(:,1)); 

CtrlVar.MinSurfAccRequiredToReactivateNodes=1e10;
iDeactivatedElements=FindElementsToDeactivate(CtrlVar,MUA1.connectivity,MUA1.coordinates,h1);
[coordinates,connectivity]=DeactivateElements(CtrlVar,iDeactivatedElements,MUA1.coordinates,MUA1.connectivity);
clear MUA1
MUA1=CreateMUA(CtrlVar,connectivity,coordinates);
figure ;  PlotFEmesh(MUA1.coordinates,MUA1.connectivity)

PlotBoundary(MUA1.Boundary,MUA1.connectivity,MUA1.coordinates,CtrlVar)

%%
% grid 2
MeshBoundaryCoordinates=[-1 -1 ; -1 1 ; 1 1 ; 1 -1];
MUA2=genmesh2d(CtrlVar,MeshBoundaryCoordinates);
CtrlVar.MeshColor='g';
CtrlVar.PlotMesh=1;
PlotFEmesh(MUA2.coordinates,MUA2.connectivity,CtrlVar)

% define a variable on the mesh1


%%

figure
PlotNodalBasedQuantities(MUA1.connectivity,MUA1.coordinates,h1);
% define a 

%%
clc
OutsideValue=-999; 
x2=[0 -1  1 0 -0.6 -0.6 -0.6 -0.6 NaN ]' ; 
y2=[0 -1 -1 1  0.0  0.4  0.6  0.7 NaN]'; 
h1=MUA1.coordinates(:,1); 
s1=h1+10;
rho1=s1*0+900;

[h2,s2,rho2]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUA1,x2,y2,OutsideValue,h1,s1,rho1);
h2'
s2'
rho2'


