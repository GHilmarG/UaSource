%% examples of how to generate FE mesh using genmesh2d 
%
% When usign Ua, the call to genmesh2d is not needed as this call is done from within Ua
%
%% example: simple square
% mesh boundary coordinates should go clockwise around the domain
% (although if this is done incorrectly, Ua will automatically correct for this anyhow.)
%
CtrlVar=Ua2D_DefaultParameters(); CtrlVar.MeshSizeMax=0.1; CtrlVar.MeshSizeMin=0.1;
MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 0];

[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)


%% example: periodic boundary conditions
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





%% example: sinusoidal bed

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

%% example% ice dome


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



%% example:  house withouth a window
CtrlVar=Ua2D_DefaultParameters(); CtrlVar.MeshSizeMax=0.1; CtrlVar.MeshSizeMin=0.1;
MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ] ;
[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

%% example: a house with a window! :-)
%
% when creating holes within a mesh, separate boundaries by NaN NaN
%
CtrlVar=Ua2D_DefaultParameters(); CtrlVar.MeshSizeMax=0.1; CtrlVar.MeshSizeMin=0.1;
MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ; ...       % Outer boundary (clockwise orientation)
               NaN NaN ;  0.5 -0.5 ; 0.5 0 ; 0.1 0 ; 0.1 -0.5 ; ...         % inner boundary (anticlockwise orientation)
               NaN NaN ; -0.1 -0.5 ; -0.1 0 ; -0.8 0 ; -0.8 -0.5 ];         % another innner boundary (anticlockwise orientation)
[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

%% a house with a tree
%
% When generating separate meshed domains, label each domain with a number.
% The label is the specified by puttin `Label NaN' ahead of the corresponding boundary
CtrlVar=Ua2D_DefaultParameters(); CtrlVar.MeshSizeMax=0.1; CtrlVar.MeshSizeMin=0.1;
MeshBoundaryCoordinates=[1 NaN ; -1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ; ...                                                     % boundary of mesh 1
                         2 NaN ; -2.0 -0.5 ; -2.0 0.5 ; -1.5 0.5 ; -1.5 -0.5 ; -1.7 -0.5 ; -1.7 -1 ; -1.8 -1 ; -1.8 -0.5 ];       % boundary of mesh 2
[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)


%% a house with windows and a tree
% When generating several separate meshed domains with holes, make sure
% to specify the outer boundary first, then the holes, and indicate to which mesh a given boundary belongs
CtrlVar=Ua2D_DefaultParameters(); CtrlVar.MeshSizeMax=0.1; CtrlVar.MeshSizeMin=0.1;
MeshBoundaryCoordinates=[1 NaN ; -1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ; ...      % outer boundary of mesh 1 (clockwise)
    1 NaN ; 0.5 -0.5 ; 0.5 0 ; 0.1 0 ; 0.1 -0.5 ; ...                              % a hole within mesh 1     (anticlockwise)
    1 NaN ; -0.1 -0.5 ; -0.1 0 ; -0.8 0 ; -0.8 -0.5 ; ...                          % a further hole within mesh 1 (anticlockwise)
    3 NaN ; -2.0 -0.5 ; -2.0 0.5 ; -1.5 0.5 ; -1.5 -0.5 ; -1.7 -0.5 ; -1.7 -1 ; -1.8 -1 ; -1.8 -0.5 ];   % another mesh (clockwise). mesh labels do not need to be in sequential order
[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

%% mesh with several holes and islands
CtrlVar=Ua2D_DefaultParameters(); CtrlVar.MeshSizeMax=0.1; CtrlVar.MeshSizeMin=0.1; CtrlVar.TriNodes=10 ;
MeshBoundaryCoordinates=...
    [1 NaN ; -1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ; ...     % outer boundary of mesh 1
    1 NaN ; 0.5 -0.5 ; 0.5 0 ; 0.1 0 ; 0.1 -0.5 ; ...         % inner boundary of mesh 1
    1 NaN ; -0.1 -0.5 ; -0.1 0 ; -0.8 0 ; -0.8 -0.5 ; ...     % another inner boundary of mesh 1
    2 NaN ; -0.7 -0.4 ; -0.7 -0.1 ; -0.2 -0.1 ; -0.2 -0.4 ; ...   % outer boundary of mesh 2
    40 NaN ; -3.0 -1.0 ; -3.0  0.5 ; -1.5 0.5 ; -1.5 -1.0 ;...    % outer boundary of mesh 40
    40 NaN ; -2.8 -0.8 ; -1.8 -0.8 ; -1.8 0.4 ; -2.8  0.4 ; ...   % inner boundary of mesh 40 
    50 NaN ; -2.6 -0.6 ; -2.6 0.3 ; -2.0 0.3 ; -2.0 -0.6  ; ...   % outer boundary of mesh 50
    50 NaN ; -2.5 -0.5 ; -2.1 -0.5 ; -2.1 0.1 ; -2.5 0.1 ];       % inner boundary of mesh 50

CtrlVar.GmeshMeshingAlgorithm=8;    % see gmsh manual


[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

% also calculate and plot normals
[nx,ny,xn,yn,Nx,Ny] = CalcEdgeAndNodalNormals(MUA.connectivity,MUA.coordinates,MUA.Boundary.Edges);
QuiverColorGHG(MUA.coordinates(MUA.Boundary.Nodes,1),MUA.coordinates(MUA.Boundary.Nodes,2),...
     Nx(MUA.Boundary.Nodes),Ny(MUA.Boundary.Nodes),[]);

%uiverColorGHG(MUA.coordinates(:,1),MUA.coordinates(:,2),Nx,Ny,[]);  % This would also work

%% example of running gmsh directly for a given input file
status=system('C:\cygwin64\home\Hilmar\GHG\ssa\Ua2D_Development\gmsh-2.8.4-Windows\gmsh.exe GmeshFileHole.geo -2 -v 1');
Gmesh=load_gmshGHG('GmeshFileHole.msh'); % the .msh is a gmsh output file. This file must have be generated previously 
TRI=Gmesh.TRIANGLES(1:Gmesh.nbTriangles,1:3);
xy=Gmesh.POS(1:Gmesh.nbNod,1:2);
figure
triplot(TRI,xy(:,1),xy(:,2)) ; axis equal

