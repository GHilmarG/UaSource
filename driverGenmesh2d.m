%% examples of how to generate FE mesh using genmesh2d 
%% example 0
CtrlVar=Ua2D_DefaultParameters();
CtrlVar.MeshSizeMax=0.1;
CtrlVar.MeshSizeMin=0.1;
MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 0];

[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

%% example for just running gmsh directly for a given input file
status=system('C:\cygwin64\home\Hilmar\GHG\ssa\Ua2D_Development\gmsh-2.8.4-Windows\gmsh.exe GmeshFileHole.geo -2 -v 1');
Gmesh=load_gmshGHG('GmeshFileHole.msh');
TRI=Gmesh.TRIANGLES(1:Gmesh.nbTriangles,1:3);
xy=Gmesh.POS(1:Gmesh.nbNod,1:2);
figure
triplot(TRI,xy(:,1),xy(:,2)) ; axis equal

%% example 1 : periodic boundary conditions
CtrlVar=Ua2D_DefaultParameters();
L=20e3 ; H=1e3;
CtrlVar.MeshSizeMax=H/10;
CtrlVar.MeshSizeMin=CtrlVar.MeshSizeMax/10;
MeshBoundaryCoordinates=[0 0 ; 0 H ; L H ; L 0];
% these lines are added to the gmsh .geo input file each time such a file is created
CtrlVar.GmshGeoFileAdditionalInputLines{1}='Physical Line(1) = {1};';  
CtrlVar.GmshGeoFileAdditionalInputLines{2}='Physical Line(2) = {2};';  
CtrlVar.GmshGeoFileAdditionalInputLines{3}='Physical Line(3) = {3};';  
CtrlVar.GmshGeoFileAdditionalInputLines{4}='Physical Line(4) = {4};';  
CtrlVar.GmshGeoFileAdditionalInputLines{5}='Physical Surface(1) = {1};';  
CtrlVar.GmshGeoFileAdditionalInputLines{6}='Periodic Line {1,2} = {3,4};';
% now everything needed for mesh generation has been defined

[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)





%% example 2
clear all
CtrlVar=Ua2D_DefaultParameters();
xL=-20e3 ; xR=20e3 ; yB=0 ; yT=2e3;
lambda=(xR-xL)/2; Ampl=1e3 ;
CtrlVar.MeshSize=0.5e3;
CtrlVar.MeshSizeMin=CtrlVar.MeshSize/5;
CtrlVar.MeshSizeMax=5*CtrlVar.MeshSize;
CtrlVar.GmeshCharacteristicLengthExtendFromBoundary=1;
CtrlVar.GmeshCharacteristicLengthFromCurvature = 1 ;

xbed=xL:CtrlVar.MeshSize:xR;
ybed=Ampl*sin(2*pi*xbed/lambda)+  yB;
xbed=xbed(:) ; ybed=ybed(:);

MeshBoundaryCoordinates=[xR yB ; xR yT ; xL yT ; xL yB ; xbed(2:end-1) ybed(2:end-1) ]; 
N=length(xbed)+2;
% these lines are added to the gmsh .geo input file each time such a file is created
CtrlVar.GmshGeoFileAdditionalInputLines{1}='Periodic Line {1} = {3};';  
CtrlVar.GmshGeoFileAdditionalInputLines{2}='Physical Line(0) = {1};';  
CtrlVar.GmshGeoFileAdditionalInputLines{3}='Physical Line(1) = {2};';  
CtrlVar.GmshGeoFileAdditionalInputLines{4}='Physical Line(2) = {3};';  
CtrlVar.GmshGeoFileAdditionalInputLines{5}=['Physical Line(3) = {4:',num2str(N),'};'];  
CtrlVar.GmshGeoFileAdditionalInputLines{6}='Physical Surface(1) = {1};';  

[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

%% example 3
% ice dome

clear all
CtrlVar=Ua2D_DefaultParameters();
L=20e3; h0=250;
xL=-L ; xR=L ; yB=0 ;  s0=3e3;

gamm=-s0/L^2;
CtrlVar.MeshSize=1e3;
CtrlVar.MeshSizeMin=0.1e3;
CtrlVar.MeshSizeMax=1e3;

xSurf=xL:CtrlVar.MeshSize:xR;
ySurf=s0+gamm*(abs(xSurf)).^2;
N=numel(xSurf);

%MeshBoundaryCoordinates=[L 0 ; -L 0 ; xSurf(:) ySurf(:)+100]; 
MeshBoundaryCoordinates=[xSurf(:) ySurf(:)]; 
% these lines are added to the gmsh .geo input file each time such a file is created
CtrlVar.GmshGeoFileAdditionalInputLines{1}=['Physical Line(1) = {',num2str(N),'};'];  
CtrlVar.GmshGeoFileAdditionalInputLines{2}='Physical Surface(1) = {1};';  

[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)



%% example:  mesh with a hole
CtrlVar=Ua2D_DefaultParameters();
CtrlVar.MeshSizeMax=0.1;
CtrlVar.MeshSizeMin=0.1;
MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ; ...
               NaN NaN ; 0.5 -0.5 ; 0.5 0 ; 0.1 0 ; 0.1 -0.5 ; ...
               NaN NaN ; -0.1 -0.5 ; -0.1 0 ; -0.8 0 ; -0.8 -0.5 ];

[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)


