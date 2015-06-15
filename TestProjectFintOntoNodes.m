CtrlVar=Ua2D_DefaultParameters();
L=20e3 ; H=20e3;
CtrlVar.MeshSizeMax=H/10;
CtrlVar.MeshSizeMin=CtrlVar.MeshSizeMax/10;
MeshBoundaryCoordinates=[0 0 ; 0 H ; L H ; L 0];
% these lines are added to the gmsh .geo input file each time such a file is created

[MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotFEmesh(MUA.coordinates,MUA.connectivity)


%%
clc
AGlen=zeros(MUA.Nnodes,1)+1 ; n=1;
x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2);
u=zeros(MUA.Nnodes,1)+sin(2*pi*x/L)+sin(2*pi*y/L);
v=zeros(MUA.Nnodes,1)+sin(2*pi*y/L)+sin(2*pi*x/L);   

[etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,u,v,AGlen,n);

[eyyNod,exxNod,exyNod]=ProjectFintOntoNodes(MUA,eyy,exx,exy);

figure; [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,u,CtrlVar)    ;
title('u')

figure; [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,v,CtrlVar)    ;
title('v')

figure; [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,exxNod,CtrlVar)    ;
title('exx')

residuals=exxNod-2*pi*cos(2*pi*x/L)/L;
figure; [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,residuals,CtrlVar)    ;
title('exx residuals')


figure; [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,eyyNod,CtrlVar)    ;
title('eyy')

residuals=eyyNod-2*pi*cos(2*pi*y/L)/L;
figure; [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,residuals,CtrlVar)    ;
title('eyy residuals')

figure; [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,exyNod,CtrlVar)    ;
title('exy')

residuals=exyNod-pi*cos(2*pi*y/L)/L-pi*cos(2*pi*x/L)/L;
figure; [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,residuals,CtrlVar)    ;
title('exy residuals')

