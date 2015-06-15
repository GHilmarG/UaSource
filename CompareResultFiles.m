
function CompareResultFiles(FileName)

if nargin==0
        [FileName1,PathName,FilterIndes]=uigetfile('*.mat');
        [FileName2,PathName,FilterIndes]=uigetfile('*.mat');
end

if isequal(FileName1,0) || isequal(FileName2,0) ; return ; end

%% set path
locdir=pwd;

indsGHG=strfind(upper(locdir),'GHG');
% 
% addpath([locdir(1:indsGHG+2),'/MatlabPackages'],'-begin')
% addpath([locdir(1:indsGHG+2),'/MatlabPackages/Mesh2d/Mesh2d v24'],'-begin')
addpath([locdir(1:indsGHG+2),'/my_mathlab_functions'],'-begin')
% addpath([locdir(1:indsGHG+2),'/ssa/FEicestream2d'],'-begin')

%%

load(FileName1,'u','v','coordinates','connectivity','nip','AGlen','n','CtrlVar','rho','rhow','h','s','S','B','b','DTxy','TRIxy','MeshBoundaryCoordinates','dhdt','C','m','n','time','CtrlVar')
u1=u ; v1=v; h1=h ; s1=s ; S1=S ; B1=B ; b1=b; dhdt1=dhdt;
load(FileName2,'u','v','coordinates','connectivity','nip','AGlen','n','CtrlVar','rho','rhow','h','s','S','B','b','DTxy','TRIxy','MeshBoundaryCoordinates','dhdt','C','m','n','time','CtrlVar')
u2=u ; v2=v; h2=h ; s2=s ; S2=S ; B2=B ; b2=b; dhdt2=dhdt;

CtrlVar=SSS2dDefaultParameters();
%[Experiment,CtrlVar,time,dt,MeshBoundaryCoordinates]=SSS2dInitialUserInput(CtrlVar);
%[DTxy,TRIxy]=TriangulationNodesIntegrationPoints(coordinates,connectivity,MeshBoundaryCoordinates);

%Smoothness=0.01 ; dx=5e3 ; dy=dx; h=SmoothFE2dMeshVariables(coordinates,h,Smoothness,dx,dy);
%CtrlVar.PlotStrains=1;
rho=rho+zeros(length(coordinates),1);

nip=1;
[etaInt,xint,yint,exx,eyy,exy,Eint,e1]=calcStrainRatesEtaInt(u1,v1,coordinates,connectivity,nip,AGlen,n,CtrlVar);
[w1,wint1]=calcVerticalSurfaceVelocity(rho,rhow,h,S1,B1,b1,u1,v1,exx,eyy,xint,yint,coordinates,connectivity,nip,CtrlVar);
[etaInt,xint,yint,exx,eyy,exy,Eint,e2]=calcStrainRatesEtaInt(u2,v2,coordinates,connectivity,nip,AGlen,n,CtrlVar);
[w2,wint2]=calcVerticalSurfaceVelocity(rho,rhow,h2,S2,B2,b2,u2,v2,exx,eyy,xint,yint,coordinates,connectivity,nip,CtrlVar);

GF = GL2d(B,S,h,rhow,rho,connectivity,CtrlVar);  % gf is only needed for plotting purposes and global remeshing
CtrlVar.PlotBackgroundImage=0;
FE2dPlots(DTxy,TRIxy,MeshBoundaryCoordinates,GF,coordinates,connectivity,b,B,S,s2-s1,h2-h1,u2-u1,v2-v1,w2-w1,dhdt2-dhdt1,C,AGlen,m,n,xint,yint,wint2-wint1,etaInt,exx,eyy,exy,e2-e1,time,CtrlVar);



