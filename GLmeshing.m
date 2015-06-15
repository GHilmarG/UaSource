% GLmeshing


load TestSave

GF = GL2d(B,S,h,rhow,rho,connectivity,CtrlVar);  % gf is only needed for plotting purposes and global remeshing



figure ; PlotFEmesh(coordinates,connectivity,CtrlVar)
hold on

CtrlVar.GLthreshold=0.5;
[GLgeo,GLinfo]=GLgeometry(connectivity,coordinates,GF,CtrlVar.GLthreshold);
plot(GLgeo(:,[3 4])',GLgeo(:,[5 6])','r')

CtrlVar.GLthreshold=0.1;
[GLgeo,GLinfo]=GLgeometry(connectivity,coordinates,GF,CtrlVar.GLthreshold);
plot(GLgeo(:,[3 4])',GLgeo(:,[5 6])','g')

%%
CtrlVar.MeshSize=10e3; CtrlVar.TriNodes=3;
[coordinates2,connectivity2]=genmesh2d(Experiment,MeshBoundaryCoordinates,CtrlVar);
figure ; PlotFEmesh(coordinates2,connectivity2,CtrlVar)

%%
xGL=GLgeo(:,7);  yGL=GLgeo(:,8); 
[yGL,ind]=sort(yGL) ; xGL=xGL(ind) ;

[xGL,yGL] = Arrange2dPos(xGL,yGL);

figure ; 
plot(xGL,yGL)

ds=100; smoothing=0.01; 
[xGL,yGL] = Smooth2dPos(xGL,yGL,ds,smoothing);

hold on
plot(xGL,yGL,'-ro')



X=[MeshBoundaryCoordinates; [xGL(:) yGL(:)]];
%X=[0 0 ; 2 0 ; 2 1 ; 0 1 ; [xGL' yGL']];
iGL=[4+(1:(numel(xGL)))];
GLedge=[iGL(1:end-1)' iGL(2:end)'];

edge=[1 iGL(1) ; iGL(1) 2 ; 2 3 ; 3 iGL(end) ; iGL(end) 4 ; 4 1 ; GLedge ];
face{1}=[1;[7:6+size(GLedge,1)]';5;6];
face{2}=[2;3;4;[6+size(GLedge,1):-1:7]'];

hdata.MeshSizeMax=10e3;
options.output=false;
%[p,t,stats] = mesh2d(X,edge) ; %,hdata,options);
[coordinates,connectivity]=genmesh2d(Experiment,X,CtrlVar,edge,face);

%[coordinates,connectivity,fnum,stats]=meshfaces(X,edge,face,hdata,options);



