%% Generating FE meshes from with Úa
%
% Several examples of how to define meshes
%
% Note that when using Úa, the call to genmesh2d is not needed as this call is done from within Úa
% Also the call CtrlVar=Ua2D_DefaultParameters() is not needed either.
%
% To run individual examples you can use the matlab option of running code sections from within editor. 
% See: 'doc run code sections'
% Just click on some part of the code using the mouse, the text will be highligted in yellow, then press ctrl ret
%
%% Example: A simple polygon 
% mesh boundary coordinates should go clockwise around the domain
% (although if this is done incorrectly, Úa will automatically correct for this anyhow.)
%
UserVar=[];
CtrlVar=Ua2D_DefaultParameters(); %
% Note; When creating this mesh using Úa, only the following 
% three lines are required in the Ua2D_InitialUserInput.m

CtrlVar.MeshGenerator='mesh2d';
CtrlVar.MeshSizeMax=0.1; 
CtrlVar.MeshSize=0.1; 
CtrlVar.MeshSizeMin=0.1;
CtrlVar.GmshMeshingAlgorithm=8;    % see gmsh manual
MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 0];
CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
% Now generate mesh (When using Úa this is done internally, no such call
% then needed).
opts.kind = 'delfront';
opts.rho2 = 1.0 ;
hfun=CtrlVar.MeshSizeMax;

tic
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates);
%[coordinates,edge,connectivity] = refine2(MeshBoundaryCoordinates,[],[],opts,hfun);
%MUA=CreateMUA(CtrlVar,connectivity,coordinates);
toc

figure ;  
subplot(1,2,1)
PlotFEmesh(MUA.coordinates,MUA.connectivity)

CtrlVar.MeshGenerator='gmsh';
tic
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates);
toc

subplot(1,2,2)
PlotFEmesh(MUA.coordinates,MUA.connectivity)

hold on 
CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
GmshPointsInLoop=PlotGmshGeometryDefinition(CtrlVar);
drawnow
%str=input('Next example? y/n [y] ? ','s');  if strcmpi(str,'n') ; return ; end
%% Example: periodic boundary conditions
CtrlVar=Ua2D_DefaultParameters();
CtrlVar.MeshGenerator='mesh2d';
UserVar=[];
% When using Úa only the following lines are needed in the input file
% Ua2D_InitialUserInput.m
L=5e3 ; H=1e3;
CtrlVar.MeshSizeMax=H/5;
CtrlVar.MeshSize=H/5;
CtrlVar.MeshSizeMin=CtrlVar.MeshSizeMax/5;
MeshBoundaryCoordinates=[0 0 ; 0 H ; L H ; L 0];
% If we want to use periodic boundary conditions then we must make sure that the
% mesh has the same periodic properties.
% These lines are added to the gmsh .geo input file each time such a file is
% created. This is a direct input to gmsh and must be on a correct from. The
% gmsh manual (http://gmsh.info/doc/texinfo/gmsh.html) gives a good description. 
% `CtrlVar.GmshGeoFileAdditionalInputLines' is a cell array with each element of
% that array being a string. This is added to the gmsh input file. This allows
% for any additional inputs to gmsh to be specified. 
CtrlVar.GmshGeoFileAdditionalInputLines{1}='Physical Line(1) = {1};';  
CtrlVar.GmshGeoFileAdditionalInputLines{2}='Physical Line(2) = {2};';  
CtrlVar.GmshGeoFileAdditionalInputLines{3}='Physical Line(3) = {3};';  
CtrlVar.GmshGeoFileAdditionalInputLines{4}='Physical Line(4) = {4};';  
CtrlVar.GmshGeoFileAdditionalInputLines{5}='Physical Surface(1) = {1};';  
CtrlVar.GmshGeoFileAdditionalInputLines{6}='Periodic Line {1,2} = {3,4};';
% Now everything needed for mesh generation has been defined. If we are running
% Úa this is all we need to do. But to see the resulting mesh we now
% generate the mesh in the exact same way as Úa would do, i.e. through a call to
% 'genmesh2d'.
UserVar=[];
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
GmshPointsInLoop=PlotGmshGeometryDefinition(CtrlVar);
drawnow

%str=input('Next example? y/n [y] ? ','s');  if strcmpi(str,'n') ; return ; end

%% Example: sinusoidal bed

CtrlVar=Ua2D_DefaultParameters();
UserVar=[];
xL=-20e3 ; xR=20e3 ; yB=0 ; yT=2e3;
lambda=(xR-xL)/2; Ampl=1e3 ;
CtrlVar.MeshSize=0.5e3;
CtrlVar.MeshSizeMin=CtrlVar.MeshSize/5;
CtrlVar.MeshSizeMax=5*CtrlVar.MeshSize;
% 
CtrlVar.GmshCharacteristicLengthExtendFromBoundary=1;
CtrlVar.GmshCharacteristicLengthFromCurvature = 1 ;

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
UserVar=[];
[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)



CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
GmshPointsInLoop=PlotGmshGeometryDefinition(CtrlVar);
drawnow

%str=input('Next example? y/n [y] ? ','s'); if strcmpi(str,'n') ; return ; end
%% Example: ice dome


CtrlVar=Ua2D_DefaultParameters();
CtrlVar.MeshGenerator='mesh2d';
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
UserVar=[];

[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)
drawnow
%str=input('Next example? y/n [y] ? ','s'); if strcmpi(str,'n') ; return ; end

%% Example:  house without a window
CtrlVar=Ua2D_DefaultParameters();
CtrlVar.GmshMeshingAlgorithm=8;    % see gmsh manual
CtrlVar.MeshSizeMax=0.1; 
CtrlVar.MeshSize=0.1; 
CtrlVar.MeshSizeMin=0.1;
MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ] ;
UserVar=[];
[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates); 
figure ; 
subplot(2,1,1)
PlotFEmesh(MUA.coordinates,MUA.connectivity)

CtrlVar.MeshGenerator='mesh2d';
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates); 
subplot(2,1,2);
PlotFEmesh(MUA.coordinates,MUA.connectivity)

drawnow
%str=input('Next example? y/n [y] ? ','s'); if strcmpi(str,'n') ; return ; end
%% Example: A house with a window! :-)
%
% when creating holes within a mesh, separate boundaries by NaN NaN
%
CtrlVar=Ua2D_DefaultParameters(); 
CtrlVar.MeshSizeMax=0.1; 
CtrlVar.MeshSize=0.1; 
CtrlVar.MeshSizeMin=0.1;


CtrlVar.MeshGenerator='mesh2d';

MeshBoundaryCoordinates=[-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ; ...       % Outer boundary (clockwise orientation)
               NaN NaN ;  0.5 -0.5 ; 0.5 0 ; 0.1 0 ; 0.1 -0.5 ; ...         % inner boundary (anticlockwise orientation)
               NaN NaN ; -0.1 -0.5 ; -0.1 0 ; -0.8 0 ; -0.8 -0.5 ];         % another innner boundary (anticlockwise orientation)
UserVar=[];           
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)
drawnow
%str=input('Next example? y/n [y] ? ','s'); if strcmpi(str,'n') ; return ; end
%% Example: A house with a tree
%
% When generating separate meshed domains, label each domain with a number.
% The label is the specified by putting `Label NaN' ahead of the corresponding boundary
CtrlVar=Ua2D_DefaultParameters();
CtrlVar.MeshGenerator='mesh2d';
CtrlVar.MeshSizeMax=0.1; 
CtrlVar.MeshSize=0.1; 
CtrlVar.MeshSizeMin=0.1; 
UserVar=[];
MeshBoundaryCoordinates=[1 NaN ; -1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ; ...                                                     % boundary of mesh 1
                         2 NaN ; -2.0 -0.5 ; -2.0 0.5 ; -1.5 0.5 ; -1.5 -0.5 ; -1.7 -0.5 ; -1.7 -1 ; -1.8 -1 ; -1.8 -0.5 ];       % boundary of mesh 2
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

drawnow
%str=input('Next example? y/n [y] ? ','s'); if strcmpi(str,'n') ; return ; end
%% Example: A house with windows and a tree :-)
% When generating several separate meshed domains with holes, make sure
% to specify the outer boundary first, then the holes, and indicate to which mesh a given boundary belongs
CtrlVar=Ua2D_DefaultParameters(); 
CtrlVar.MeshSizeMax=0.1; 
CtrlVar.MeshSizeMin=0.1;
CtrlVar.MeshSize=0.1;
CtrlVar.MeshGenerator='mesh2d';
UserVar=[];
MeshBoundaryCoordinates=[1 NaN ; -1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ; ...      % outer boundary of mesh 1 (clockwise)
    1 NaN ; 0.5 -0.5 ; 0.5 0 ; 0.1 0 ; 0.1 -0.5 ; ...                              % a hole within mesh 1     (anticlockwise)
    1 NaN ; -0.1 -0.5 ; -0.1 0 ; -0.8 0 ; -0.8 -0.5 ; ...                          % a further hole within mesh 1 (anticlockwise)
    3 NaN ; -2.0 -0.5 ; -2.0 0.5 ; -1.5 0.5 ; -1.5 -0.5 ; -1.7 -0.5 ; -1.7 -1 ; -1.8 -1 ; -1.8 -0.5 ];   % another mesh (clockwise). mesh labels do not need to be in sequential order
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)
drawnow
%str=input('Next example? y/n [y] ? ','s'); if strcmpi(str,'n') ; return ; end
%% Example: Mesh with several holes and islands
CtrlVar=Ua2D_DefaultParameters(); 
CtrlVar.MeshSizeMax=0.1; 
CtrlVar.MeshSizeMin=0.1; 
CtrlVar.MeshSize=0.1; 
CtrlVar.TriNodes=10 ;
CtrlVar.MeshGenerator='mesh2d';
MeshBoundaryCoordinates=...
    [1 NaN ; -1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ; ...     % outer boundary of mesh 1
    1 NaN ; 0.5 -0.5 ; 0.5 0 ; 0.1 0 ; 0.1 -0.5 ; ...         % inner boundary of mesh 1
    1 NaN ; -0.1 -0.5 ; -0.1 0 ; -0.8 0 ; -0.8 -0.5 ; ...     % another inner boundary of mesh 1
    2 NaN ; -0.7 -0.4 ; -0.7 -0.1 ; -0.2 -0.1 ; -0.2 -0.4 ; ...   % outer boundary of mesh 2
    40 NaN ; -3.0 -1.0 ; -3.0  0.5 ; -1.5 0.5 ; -1.5 -1.0 ;...    % outer boundary of mesh 40
    40 NaN ; -2.8 -0.8 ; -1.8 -0.8 ; -1.8 0.4 ; -2.8  0.4 ; ...   % inner boundary of mesh 40 
    50 NaN ; -2.6 -0.6 ; -2.6 0.3 ; -2.0 0.3 ; -2.0 -0.6  ; ...   % outer boundary of mesh 50
    50 NaN ; -2.5 -0.5 ; -2.1 -0.5 ; -2.1 0.1 ; -2.5 0.1 ];       % inner boundary of mesh 50

CtrlVar.GmshMeshingAlgorithm=8;    % see gmsh manual

UserVar=[];
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

% also calculate and plot normals
[nx,ny,xn,yn,Nx,Ny] = CalcEdgeAndNodalNormals(MUA.connectivity,MUA.coordinates,MUA.Boundary.Edges);
QuiverColorGHG(MUA.coordinates(MUA.Boundary.Nodes,1),MUA.coordinates(MUA.Boundary.Nodes,2),...
     Nx(MUA.Boundary.Nodes),Ny(MUA.Boundary.Nodes),[]);
 
CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
GmshPointsInLoop=PlotGmshGeometryDefinition(CtrlVar);
 
drawnow

 %str=input('Next example? y/n [y] ? ','s'); if strcmpi(str,'n') ; return ; end

 
%% Example: an elliptical hole/crack
%
% when creating holes within a mesh, separate boundaries by NaN NaN
%
CtrlVar=Ua2D_DefaultParameters(); 
CtrlVar.MeshGenerator='mesh2d';
%CtrlVar.MeshGenerator='gmsh';
CtrlVar.MeshSizeMax=1e3; 
CtrlVar.MeshSize=1e3; 
CtrlVar.MeshSizeMin=0.1e3;

a=0.1e3; % horizontal radius
b=5e3; % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
t=linspace(-pi,pi,200);  t(end)=[];
xe=x0+a*cos(t);
ye=y0+b*sin(t);

CtrlVar.GmshCharacteristicLengthFromCurvature = 1 ;
CtrlVar.GmshCharacteristicLengthExtendFromBoundary=1;


MeshBoundaryCoordinates=...
    [-10e3 -10e3 ; ...
    -10e3 10e3 ; ...
    10e3 10e3 ; ...
    10e3 -10e3 ; ...
    NaN  NaN ; ...
     xe(:) ye(:)] ;

UserVar=[];           
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates); 
figure ;  
PlotFEmesh(MUA.coordinates,MUA.connectivity)
CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
GmshPointsInLoop=PlotGmshGeometryDefinition(CtrlVar);
drawnow
%% Example: Constrainted meshing with two joined meshes sharing the same boundary
%
% Since here the same line belongs to more then one mesh we need to specify
% the 'plane surface' (gmsh terminology) separatly. This is done in CtrlVar.GmshPlaneSurface

CtrlVar=Ua2D_DefaultParameters(); 
CtrlVar.MeshSizeMax=0.1; 
CtrlVar.MeshSize=0.1; 
CtrlVar.MeshSizeMin=0.1;
CtrlVar.MeshGenerator='mesh2d';
CtrlVar.MeshGenerator='gmsh';

UserVar=[];
MeshBoundaryCoordinates=[1 NaN ;-1 -1 ; -1 0 ; 0 1 ; 1 0 ; 1 -1 ; 0 -1 ; ...       % Outer boundary (clockwise orientation)
                         1 NaN ;  0.5 -0.5 ; 0.5 0 ; 0.1 0 ; 0.1 -0.5 ; ...         % inner boundary (anticlockwise orientation)
                         1 NaN ; -0.1 -0.5 ; 0 0.3 ; -0.95 0 ; -0.8 -0.9 ];         % another innner boundary (anticlockwise orientation)

CtrlVar.GmshPlaneSurface{1}=[1 2 3]; % mesh 1 is defined by boundaries 1, 2 and 3
CtrlVar.GmshPlaneSurface{2}=[-2];  % mesh 2 is defined by boundary 2, make sure to use a negative sign to indicate a reverse ordering
CtrlVar.GmshPlaneSurface{3}=[-3];  % mesh 2 is defined by boundary 2, make sure to use a negative sign to indicate a reverse ordering
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

plot(MeshBoundaryCoordinates(:,1),MeshBoundaryCoordinates(:,2),'-or','LineWidth',2)

drawnow
%str=input('Next example? y/n [y] ? ','s'); if strcmpi(str,'n') ; return ; end

%% Example: two separate meshes in contact


CtrlVar=Ua2D_DefaultParameters(); CtrlVar.MeshSizeMax=0.1; CtrlVar.MeshSizeMin=0.1;
UserVar=[];
MeshBoundaryCoordinates=[1 NaN ;  0 0  ; 0 1 ; 1 1 ; 1 0 ; ...      
                         2 NaN ; -1 0  ; -1 1.001 ; 0 1.001 ; 0 0 ];

[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates); figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

plot(MeshBoundaryCoordinates(:,1),MeshBoundaryCoordinates(:,2),'-or','LineWidth',2)
drawnow
%str=input('Next example? y/n [y] ? ','s'); if strcmpi(str,'n') ; return ; end

%% Example: Internal boundaries


CtrlVar=Ua2D_DefaultParameters(); CtrlVar.MeshSizeMax=0.1; CtrlVar.MeshSizeMin=0.1; CtrlVar.GmshInputFormat=1;
UserVar=[];
MeshBoundaryCoordinates=[1 NaN ;  0 0  ; 0 0.25 ; 0.25 0.25 ; 0.25 0.75 ; 0 0.75 ; 0 1 ; 1 1 ; 1 0;   ...      
                         2 NaN ;  0 0.25 ; 0. 0.75 ; 0.25 0.75 ; 0.25 0.25 ];

[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates); 
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

plot(MeshBoundaryCoordinates(:,1),MeshBoundaryCoordinates(:,2),'-or','LineWidth',2)
drawnow
%str=input('Next example? y/n [y] ? ','s');  if strcmpi(str,'n') ; return ; end

%% Example: Two subdomains sharing a common boundary.
%
% This example uses GmshInputFormat 2.
% 
% Gmsh input format 2 is quite close to the gmsh input format itself and more
% flexible than the default input format 1, but also more difficult to use.
%
% The gmsh input format 2 does not use `MeshBoundaryCoordinates' at all.
%
% Using gmsh input format 2 allows constraint meshing. The computational domain can be
% split up into several subdomains. each with their own boundaries. One can create
% internal boundaries within a domain, and subdomains sharing some common boundaries.
%
% All inputs are defined as fields to the CtrlVar.
%
% The basic idea is to define points, then to specify lines in terms of those
% points, loops in terms of those lines, and finally plane surfaces in terms of
% the loops.
%
%
CtrlVar=Ua2D_DefaultParameters(); 
UserVar=[];
CtrlVar.GmshInputFormat=2; % using input format 2
CtrlVar.MeshSizeMax=0.1; 
CtrlVar.MeshSizeMin=0.1;
CtrlVar.MeshSize=0.1;

CtrlVar.Gmsh.Points=[0 0 ;...    % now define a few gmsh points
                     0 0.5 ; ... % although not directly labeled
                     0 1 ;   ... % the first point in point 1, etc. 
                     1 1 ;  ...
                     1 0 ; ... 
                     -1 0 ; ...   % this is point nr 6
                     -1 0.5];     % and this point nr 7

CtrlVar.Gmsh.Lines{1}=[1 , 2];  % define gmsh lines using the point labels
CtrlVar.Gmsh.Lines{2}=[2 , 3];
CtrlVar.Gmsh.Lines{3}=[ 3 ; 4 ; 5 ; 1 ];
CtrlVar.Gmsh.Lines{4}=[ 1 ;  6 ; 7 ; 2 ];

CtrlVar.Gmsh.Loops{1}=[1 ; 2 ; 3 ]; % now define loops
CtrlVar.Gmsh.Loops{2}=[4 ; -1];     % all loops must be closed!
                                    % (negative sign reverses the line direction.)

CtrlVar.Gmsh.PlaneSurfaces{1} = [1]; % finally define the surfaces 
CtrlVar.Gmsh.PlaneSurfaces{2} = [2];

MeshBoundaryCoordinates=[]; % When CtrlVar.GmshInputFormat=2, 'MeshBoundaryCoordinates' are not used
[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates); 
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)
hold on ; PlotGmshGeometryDefinition(CtrlVar);
drawnow
%str=input('Next example? y/n [y] ? ','s');  if strcmpi(str,'n') ; return ; end
%% Example: Several different subdomains
%
% This example uses gmsh input format 2
%
CtrlVar=Ua2D_DefaultParameters(); 
UserVar=[];
MeshBoundaryCoordinates=[];
CtrlVar.GmshInputFormat=2;
CtrlVar.MeshSizeMax=0.1; CtrlVar.MeshSizeMin=0.1; CtrlVar.TriNodes=3;

CtrlVar.Gmsh.Points=[0 0 ;...  
                     0 0.5 ; ...
                     0 1 ;   ...
                     1 1 ;  ...
                     1 0 ; ... 
                     -1 0 ; ...
                     -1 0.5 ; ...
                     0.25 0.25 ; ...
                     0.35 0.75 ; ...
                     0.65 0.75 ;
                     0.80 0.25 ;
                     -0.5 1 ; ...
                     -1 2 ; ...
                     1 2 ; ...
                     ];

CtrlVar.Gmsh.Lines{1}=[1 ; 2];
CtrlVar.Gmsh.Lines{2}=[2 ; 3];
CtrlVar.Gmsh.Lines{3}=[ 3 ; 4 ; 5 ; 1 ];
CtrlVar.Gmsh.Lines{4}=[ 1 ;  6 ; 7 ; 2 ];
CtrlVar.Gmsh.Lines{5}=[ 8 ; 9 ; 10 ; 11 ; 8];
CtrlVar.Gmsh.Lines{6}=[ 12 ; 13 ; 14 ; 12];

CtrlVar.Gmsh.Loops{1}=[1 ; 2 ; 3 ];
CtrlVar.Gmsh.Loops{2}=[4 ; -1];
CtrlVar.Gmsh.Loops{3}=[-5]; % reverse orientation 
CtrlVar.Gmsh.Loops{4}=[6];
CtrlVar.Gmsh.Loops{5}=[5];  %

CtrlVar.Gmsh.PlaneSurfaces{1} = [1 ; 3 ]; % this surface has a hole
CtrlVar.Gmsh.PlaneSurfaces{2} = [2];
CtrlVar.Gmsh.PlaneSurfaces{3} = [4];
CtrlVar.Gmsh.PlaneSurfaces{4} = [5];


[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates); 

% plotting mesh
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)
hold on
PlotGmshGeometryDefinition(CtrlVar);

drawnow
%str=input('Next example? y/n [y] ? ','s');  if strcmpi(str,'n') ; return ; end

%% Example: Meshing the Brunt Ice Shelf
%

load BruntMeshBoundaryCoordinates.mat
UserVar=[];
I=find(isnan(MeshBoundaryCoordinates(:,1)) | isnan(MeshBoundaryCoordinates(:,2)));

i1a=I(1)+1 ; i1b=I(2)-1;  
i2a=I(2)+1 ; i2b=size(MeshBoundaryCoordinates,1);

CtrlVar=Ua2D_DefaultParameters(); 
CtrlVar.GmshInputFormat=1;
CtrlVar.MeshSizeMax=5e3; 
CtrlVar.MeshSizeMin=1e3; 
CtrlVar.MeshSize=5e3; 
CtrlVar.TriNodes=3;
CtrlVar.MeshGenerator='mesh2d';
%CtrlVar.MeshGenerator='gmsh';

CtrlVar.Gmsh.Points=MeshBoundaryCoordinates(i1a:i1b,:);

CtrlVar.Gmsh.Loops{1}=[1];
CtrlVar.Gmsh.PlaneSurfaces{1} = [1];
CtrlVar.Gmsh.Lines{1}=[(1:i1b-i1a+1)';1];
[UserVar,MUA]=genmesh2d2(UserVar,CtrlVar,MeshBoundaryCoordinates); 


% plotting mesh
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)
drawnow
%str=input('Next example? y/n [y] ? ','s');  if strcmpi(str,'n') ; return ; end
%%
% Example: Meshing the Brunt Ice Shelf with an internal boundary (gmsh input format 2)
load BruntMeshBoundaryCoordinates.mat


CtrlVar=Ua2D_DefaultParameters(); 
CtrlVar.GmshInputFormat=2;
CtrlVar.MeshSizeMax=5e3; 
CtrlVar.MeshSizeMin=1e3; 
CtrlVar.TriNodes=3;




CtrlVar.Gmsh.Points=[MeshBoundaryCoordinates(i1a:i1b,:); MeshBoundaryCoordinates(i2a:i2b,:) ];
l1a=1 ; l1b=l1a+i1b-i1a;
l2a=l1b+1 ; l2b=l2a+i2b-i2a;
CtrlVar.Gmsh.Lines{1}=[(l1a:l1b)';l1a];
CtrlVar.Gmsh.Lines{2}=[(l2a:l2b)';l2a];
CtrlVar.Gmsh.Loops{1}=[1];
CtrlVar.Gmsh.Loops{2}=[2];
CtrlVar.Gmsh.PlaneSurfaces{1} = [1,-2];
CtrlVar.Gmsh.PlaneSurfaces{2} = [2];

UserVar=[];

[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates); 
CtrlVar.Gmsh.Lines{1}=[(1:i1b-i1a+1)';1];

% plotting mesh
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

hold on ; PlotGmshGeometryDefinition(CtrlVar);
drawnow
%str=input('Next example? y/n [y] ? ','s');  if strcmpi(str,'n') ; return ; end
%% Example: Meshing the Brunt Ice Shelf with an internal boundary and a `fill-in'

load BruntMeshBoundaryCoordinates.mat

I=find(isnan(MeshBoundaryCoordinates(:,1)) | isnan(MeshBoundaryCoordinates(:,2)));

i1a=I(1)+1 ; i1b=I(2)-1;  
i2a=I(2)+1 ; i2b=size(MeshBoundaryCoordinates,1);

UserVar=[];
CtrlVar=Ua2D_DefaultParameters(); 
CtrlVar.GmshInputFormat=2;
CtrlVar.MeshSizeMax=5e3; 
CtrlVar.MeshSizeMin=1e3; 
CtrlVar.TriNodes=3;
CtrlVar.Gmsh.Points=[MeshBoundaryCoordinates(i1a:i1b,:); MeshBoundaryCoordinates(i2a:i2b,:) ];
l1a=1 ; l1b=l1a+i1b-i1a;
l2a=l1b+1 ; l2b=l2a+i2b-i2a;
CtrlVar.Gmsh.Lines{1}=[(l1a:l1b)';l1a];
CtrlVar.Gmsh.Lines{2}=[(l2a:l2b)';l2a];


Box=[-7.02e+05  -6.4314e+05   1.4651e+06   1.505e+06];
I=CtrlVar.Gmsh.Points(:,1) > Box(1) & CtrlVar.Gmsh.Points(:,1) < Box(2) ...
    & CtrlVar.Gmsh.Points(:,2) > Box(3) & CtrlVar.Gmsh.Points(:,2) < Box(4) ;

% now must line 1 into two lines, with one convering the BIS-SW chasm outlines only

x=CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{1},1);
y=CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{1},2);
K=x > Box(1) & x < Box(2) ...
    & y > Box(3) & y < Box(4) ;

temp=CtrlVar.Gmsh.Lines{1};
temp(K)=NaN; temp(end)=[];
J=find(isnan(temp));
tt=[temp(J(end)+1:end,end);temp(1:J(1)-1)];
CtrlVar.Gmsh.Lines{3}=tt;
CtrlVar.Gmsh.Lines{4}=find(K);
CtrlVar.Gmsh.Lines{4}=[CtrlVar.Gmsh.Lines{3}(end); CtrlVar.Gmsh.Lines{4};CtrlVar.Gmsh.Lines{3}(1)];

CtrlVar.Gmsh.Lines{5}=[CtrlVar.Gmsh.Lines{3}(end);CtrlVar.Gmsh.Lines{3}(1)];
CtrlVar.Gmsh.Lines{6}=[CtrlVar.Gmsh.Lines{4}(end);CtrlVar.Gmsh.Lines{4}(1)];

figure 
plot(CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{1},1),CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{1},2),'.-')
hold on
plot(CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{2},1),CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{2},2),'^-')
plot(CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{3},1),CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{3},2),'+-')
hold on
plot(CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{4},1),CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{4},2),'O-')
plot(CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{5},1),CtrlVar.Gmsh.Points(CtrlVar.Gmsh.Lines{5},2),'*-')
legend('1','2','3','4','5')


CtrlVar.Gmsh.Loops{1}=[3,4];
CtrlVar.Gmsh.Loops{2}=[2];
CtrlVar.Gmsh.Loops{3}=[-4,5];

CtrlVar.Gmsh.PlaneSurfaces{1} = [1,-2];
CtrlVar.Gmsh.PlaneSurfaces{2} = [2];
CtrlVar.Gmsh.PlaneSurfaces{3} = [3];


[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates); 
CtrlVar.Gmsh.Lines{1}=[(1:i1b-i1a+1)';1];

% plotting mesh
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)

hold on ; PlotGmshGeometryDefinition(CtrlVar);
drawnow
%str=input('Next example? y/n [y] ? ','s');  if strcmpi(str,'n') ; return ; end
%% Example of running gmsh directly for a given input file

status=system([getenv('GmshHomeDirectory'),'\gmsh.exe GmshFile.geo -2 -v 5']);
Gmsh=load_gmshGHG('GmshFile.msh'); % the .msh is a gmsh output file. This file must have be generated previously 
TRI=Gmsh.TRIANGLES(1:Gmsh.nbTriangles,1:3);
xy=Gmsh.POS(1:Gmsh.nbNod,1:2);
figure
triplot(TRI,xy(:,1),xy(:,2)) ; axis equal
drawnow
%hold on
%plot(MeshBoundaryCoordinates(:,1),MeshBoundaryCoordinates(:,2),'-or','LineWidth',2)

%%