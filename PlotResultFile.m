
function PlotResultFile(FileName)

if nargin==0
        [FileName,PathName,FilterIndes]=uigetfile('*.mat');
end

if isequal(FileName,0) ; return ; end



load(FileName,'as','ab','u','v','coordinates','connectivity','nip','AGlen','n','CtrlVar','rho','rhow','g','h','s','S','B','b','DTxy','TRIxy','MeshBoundaryCoordinates','dhdt','C','m','n','time','GFstart','CtrlVar')
CtrlVar=SSS2dDefaultParameters();
%%
[a,as,ab]=DefineMassBalance(' ',time,s,b,S,B,coordinates,connectivity,TRIxy,rhow,rho,CtrlVar);

%[Experiment,CtrlVar,time,dt,MeshBoundaryCoordinates]=SSS2dInitialUserInput(CtrlVar);
%[DTxy,TRIxy]=TriangulationNodesIntegrationPoints(coordinates,connectivity,MeshBoundaryCoordinates);

%Smoothness=0.01 ; dx=5e3 ; dy=dx; h=SmoothFE2dMeshVariables(coordinates,h,Smoothness,dx,dy);
%CtrlVar.PlotStrains=1;
nip=1;
rho=rho+zeros(length(coordinates),1);

[etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
[w,wint]=calcVerticalSurfaceVelocity(rho,rhow,h,S,B,b,u,v,as,ab,exx,eyy,xint,yint,coordinates,connectivity,nip,CtrlVar);
GF = GL2d(B,S,h,rhow,rho,connectivity,CtrlVar);  % gf is only needed for plotting purposes and global remeshing
CtrlVar.PlotBackgroundImage=0;


if exist('GFstart','var')==0; GFstart=GF ; end
if numel(GF.node)==numel(GFstart.node) ;  dGFdt=(GF.node-GFstart.node); else dGFdt=u*0 ; end


FE2dPlots(CtrlVar,DTxy,TRIxy,MeshBoundaryCoordinates,GF,dGFdt,coordinates,connectivity,...
    b,B,S,s,h,u,v,w,dhdt,C,AGlen,m,n,xint,yint,wint,etaInt,exx,eyy,exy,e,time,...
    rho,rhow,as+ab,as,ab,txx,tyy,txy);


 





