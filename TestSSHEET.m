
% %% Mehs : square
% 
% CtrlVar=Ua2D_DefaultParameters();
% CtrlVar.fidlog=1;
% CtrlVar.MeshSizeMax=1e3;
% CtrlVar.MeshSizeMin=1e3;
% MeshBoundaryCoordinates=[-50e3 -50e3 ; 50e3 -50e3 ; 50e3 50e3 ; -50e3 50e3];
% CtrlVar.GmshMeshingAlgorithm=8;  % see gmsh manual 
% [coordinates,connectivity]=genmesh2d(' ',MeshBoundaryCoordinates,CtrlVar);
% 
% figure ; PlotFEmesh(coordinates,connectivity,CtrlVar)
% 
% save UAPmesh coordinates connectivity

%%
CtrlVar=Ua2D_DefaultParameters();
CtrlVar.PlotXYscale=1000;
CtrlVar.ParallelAssembly=1;
load UAPmesh coordinates connectivity
nip=6 ;
[MeshDeriv,MeshDetJ]=CalcMeshDerivatives(CtrlVar,connectivity,coordinates,nip);
Nnodes=size(coordinates,1) ; [Nele,nod]=size(connectivity);
x=coordinates(:,1) ; y=coordinates(:,2);
TestCase='gauss';

switch TestCase
    case 'uip'
        
        alpha=0.1;
        B=zeros(Nnodes,1)-tan(alpha)*coordinates(:,1);
        S=zeros(Nnodes,1)-10000;
        a=zeros(Nnodes,1) ;
        b=B;
        s=b+1000;
    case 'gauss'
        alpha=0.01;
        B=zeros(Nnodes,1)-tan(alpha)*coordinates(:,1);
        S=zeros(Nnodes,1)-10000;
        a=zeros(Nnodes,1) ;
        b=B;
        ampl=100 ; sx=10e3 ; sy=10e3;
        s=b+ 10+ ampl*exp(-(x.*x)/sx^2-(y.*y)/sy^2);
        
    case 'acc'
        alpha=0.0;
        B=zeros(Nnodes,1);
        S=zeros(Nnodes,1)-10000;
        a=zeros(Nnodes,1)+10 ;
        b=B;
        s=b+10;
        
end

rho=900+zeros(Nnodes,1);
rhow=1000;
g=9.81/1000;
secinyear=365.24*24*60*60;
AGlen=zeros(Nnodes,1)+1e-12*1000*secinyear; n=1;
CtrlVar.AGlenisElementBased=0 ;


CtrlVar.theta=0.5; 
CtrlVar.fidlog=1;



s0=s ;  b0=b ;  a1=a ; a0=a;
h0=s0-b0; 


hstart=h0;


time=0 ;
dt=10;
nTime=1;


Lh=[] ; Lhrhs=[] ; lambdah=[];

as0=a0 ; ab0=a0*0 ; as1=a1 ; ab1=a1*0;
h1=h0; b1=B; s1=b1+h1;



CtrlVar.InfoLevelNonLinIt=1;
CtrlVar.InfoLevelCPU=0;

MUA.coordinates=coordinates;
MUA.connectivity=connectivity;
MUA.nip=nip;
MUA.Deriv=MeshDeriv;
MUA.DetJ=MeshDetJ;
MUA.Nnodes=Nnodes;
MUA.Nele=Nele;
MUA.nod=nod;

[u,v]=uvSSHEET(CtrlVar,MUA,AGlen,n,rho,g,s1,h1);
CtrlVar.RelativeVelArrowSize=1e-1;
figure ; [cbar,uvPlotScale]=QuiverColorGHG(coordinates(:,1),coordinates(:,2),u,v,CtrlVar);

for I=1:nTime
    fprintf(' \n \n --------------   Time Step %-i \t time=%-g \n ',I,time)
[u1,v1,h1,s1,lambdah1,nlInfo]=SSHEET_TransientImplicit(CtrlVar,MUA,dt,s0,s1,b0,b1,as0,ab0,as1,ab1,Lh,Lhrhs,lambdah,AGlen,n,rho,g);

time=time+dt;


fprintf('\n')
figure(1) ;
subplot(2,2,1) ; PlotNodalBasedQuantities(connectivity,coordinates,hstart,CtrlVar);    title(sprintf('h at t=%-g ',0))
subplot(2,2,2) ; PlotNodalBasedQuantities(connectivity,coordinates,h1,CtrlVar);        title(sprintf('h at t=%-g ',time))
subplot(2,2,3) ; PlotNodalBasedQuantities(connectivity,coordinates,h1-h0,CtrlVar);     title(sprintf(' dh at t=%-g ',time))
subplot(2,2,4) ; PlotNodalBasedQuantities(connectivity,coordinates,h1-hstart,CtrlVar); title(sprintf(' h(t=%-g)-h(t=0) ',time))
dh=(h1-h0)/dt;
h0=h1 ; s0=s1 ;
h1=dh*dt+h0; s1=b1+h1;

end