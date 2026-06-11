%%
close all ; clear all

% call gmsh, and gmsh meshes all surfacces (-2 flag)
% and writes out a msh file
! gmsh-2.6.1-Windows/gmsh.exe t2.geo -2 

% loading msh file, returning a Gmsh structure
Gmsh=load_gmshGHG('t2.msh');

TRI=Gmsh.TRIANGLES(1:Gmsh.nbTriangles,1:3);
xy=Gmsh.POS(1:Gmsh.nbNod,1:2);
triplot(TRI,xy(:,1),xy(:,2)) ; axis equal

%% define a background scalar field an remesh area 
xTRI=reshape(xy(TRI,1)',[],3);
yTRI=reshape(xy(TRI,2)',[],3);
xyTRI=[xTRI(:,1) yTRI(:,1) xTRI(:,2) yTRI(:,2) xTRI(:,3) yTRI(:,3) ];
TargetSize=zeros(Gmsh.nbTriangles,3)+0.1;
I=xTRI>0.4 &  xTRI<0.6 ; 
TargetSize(I)=0.01;
FileName='bgmsh2.pos';
%io = CreateGmshBackgroundScalarMesh(xTRI,yTRI,TargetSize,FileName);
io = CreateGmshBackgroundScalarMesh(xy,TRI,TargetSize,FileName);

! gmsh-2.6.1-Windows/gmsh.exe t2.geo -bgm bgmsh2.pos -o t2bgm.msh -2  


Gmsh=load_gmshGHG('t2bgm.msh');

TRI=Gmsh.TRIANGLES(1:Gmsh.nbTriangles,1:3);
xy=Gmsh.POS(1:Gmsh.nbNod,1:2);
figure
triplot(TRI,xy(:,1),xy(:,2)) ; axis equal
