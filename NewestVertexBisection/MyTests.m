
%% Test 1
close all
mesh = genMesh([1,2,3;3,4,1], [0,0;1,0;1,1;0,1]) ; %, [1,2],[2,3;3,4;4,1]);
mesh = genBisectionMesh(mesh);
CtrlVar.PlotLabels=1;

mesh = bisectionRefine2D(mesh, 'all');
figure ; PlotFEmesh(mesh.coordinates,mesh.elements,CtrlVar) ;  title('mesh1')


%%
mesh2 = bisectionRefine(mesh, [1 ; 2 ;3 ]); 
figure ; PlotFEmesh(mesh2.coordinates,mesh2.elements,CtrlVar) ;   title('mesh2')

mesh3 = bisectionRefine(mesh2, [3 ; 4 ]); 
figure ; PlotFEmesh(mesh3.coordinates,mesh3.elements,CtrlVar) ; title('mesh3')



Mesh2 = bisectionCoarsen(mesh3, [4 5 8 9]); figure ; PlotFEmesh(Mesh2.coordinates,Mesh2.elements,CtrlVar) ; title('Mesh2')
Mesh1 = bisectionCoarsen(Mesh2, 'all'); figure ; PlotFEmesh(Mesh1.coordinates,Mesh1.elements,CtrlVar) ; title('Mesh1')


%% Test 2
close all
CtrlVar.PlotLabels=1;
mesh = genMesh([1,2,3; 1,3,4; 4 3 5 ; 5 6 4], [0,0;1,0;1,1;0,1; 1 2 ; 0 2]) ; %, [1,2],[2,3;3,4;4,1]);
mesh = genMesh([3,1,2 ; 4,1,3], [0,0;1,0;1,1;0,1; 1 2 ; 0 2],[1,2],[2,3;3,4;4,1]);
%mesh = genMesh([3,1,2 ; 1,3,4], [0,0;1,0;1,1;0,1; 1 2 ; 0 2],[1,2],[2,3;3,4;4,1]);

mesh=SelectRefinementEdge(mesh);

mesh.bd=[];
%mesh = genBisectionMesh(mesh,'forceRefine');
mesh = genBisectionMesh(mesh);
figure ; PlotFEmesh(mesh.coordinates,mesh.elements,CtrlVar) ;  title('mesh1')

CtrlVar.PlotLabels=1;

mesh = bisectionRefine2D(mesh, 'all');
figure ; PlotFEmesh(mesh.coordinates,mesh.elements,CtrlVar) ;  title('mesh2')


%% Test
close all
CtrlVar.PlotLabels=1;

mesh = genMesh([1,2,7 ; 2,3,7; 3 ,4 , 7 ; 4 , 5 , 7 ; 5 , 6 , 7 ; 6, 1 , 7], [-0.5,-1;0.5,-1; 1,0 ; 0.5 ,1; -0.5 1 ; -1 0 ; 0 0 ]);
mesh.bd=[]; mesh = genBisectionMesh(mesh);
figure ; PlotFEmesh(mesh.coordinates,mesh.elements,CtrlVar) ;  title('mesh1')


mesh=SelectRefinementEdge(mesh);
mesh = bisectionRefine2D(mesh, 'all');
figure ; PlotFEmesh(mesh.coordinates,mesh.elements,CtrlVar) ;  title('mesh2')

%% Test Ua
clear all
close all
load TestSave
[MUA.coordinates,MUA.connectivity]=ChangeElementType(MUA.coordinates,MUA.connectivity,3);


CtrlVar.PlotLabels=0;
%connectivity=FlipElements(MUA.connectivity);

mesh = genMesh(MUA.connectivity, MUA.coordinates);
mesh.bd=[];
mesh = genBisectionMesh(mesh);
mesh=SelectRefinementEdge(mesh);


figure ; PlotFEmesh(mesh.coordinates,mesh.elements,CtrlVar) ;  title('mesh forceRefine')


meshAll = bisectionRefine2D(mesh, 'all');
figure ; PlotFEmesh(meshAll.coordinates,meshAll.elements,CtrlVar) ;  title('mesh All')


mesh1 = bisectionRefine2D(mesh, [100:403]);
figure ; PlotFEmesh(mesh1.coordinates,mesh1.elements,CtrlVar) ;  title('mesh 1 first ')

mesh1 = bisectionRefine2D(mesh1, [200:400]);
figure ; PlotFEmesh(mesh1.coordinates,mesh1.elements,CtrlVar) ;  title('mesh 1 second')

mesh1 = bisectionRefine2D(mesh1, [200:400]);
figure ; PlotFEmesh(mesh1.coordinates,mesh1.elements,CtrlVar) ;  title('mesh 1 third')

mesh1=mesh; 
for I=1:2
    mesh1 = bisectionRefine2D(mesh1,ElementsToBeRefined);
    figure(I+100) ; PlotFEmesh(mesh1.coordinates,mesh1.elements,CtrlVar) ;  title('mesh refinement')
end

mesh1 = bisectionCoarsen(mesh1,'all') ;  %markedElements, varargin)
figure(1000) ; PlotFEmesh(mesh1.coordinates,mesh1.elements,CtrlVar) ;  title('mesh coarsen')






























