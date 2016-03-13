%%
close all ; clear all

% call gmsh, and gmsh meshes all surfacces (-2 flag)
% and writes out a msh file
! gmsh-2.6.1-Windows/gmsh.exe t1.geo -2 

Gmesh=load_gmshGHG('t1.msh');

TRI=Gmesh.TRIANGLES(1:Gmesh.nbTriangles,1:3);
xy=Gmesh.POS(1:Gmesh.nbNod,1:2);
triplot(TRI,xy(:,1),xy(:,2)) ; axis equal

%%

status=system('C:\cygwin64\home\Hilmar\ghg\Ua\gmsh-2.11.0-Windows\gmsh.exe GmeshFile.geo -2 -v 1');
Gmesh=load_gmshGHG('GmeshFile.msh');
TRI=Gmesh.TRIANGLES(1:Gmesh.nbTriangles,1:3);
xy=Gmesh.POS(1:Gmesh.nbNod,1:2);

figure
triplot(TRI,xy(:,1),xy(:,2)) ; axis equal