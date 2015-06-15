


function EleSizeIndicator=DesiredEleSizesBasedOnMeasVelocity(CtrlVar,MUA,s,b,S,B,rho,rhow,AGlen,n,GF)


[uMeas,vMeas]=EricVelocities(CtrlVar,MUA.coordinates);
e=EricStrainRates(CtrlVar,MUA,AGlen,n,uMeas,vMeas);
e=ProjectFintOntoNodes(MUA,e);
e(e<1e-3)=1e-3 ; e(e>0.025)=0.025 ;



tri=MUA.connectivity;
figure
PlotNodalBasedQuantities(tri,MUA.coordinates,e,CtrlVar);
title(' e based on measurements ')


EleSizeIndicator=real(1./(e.^CtrlVar.hpower)) ;

%CtrlVar.MeshSizeMin=2*CtrlVar.MeshSizeMin;  % Not going all the way down based on this error estimate

EleSizeIndicator=CtrlVar.MeshSizeMin+...
    (CtrlVar.MeshSizeMax-CtrlVar.MeshSizeMin)*(EleSizeIndicator-min(EleSizeIndicator))/(max(EleSizeIndicator)-min(EleSizeIndicator));

speed=sqrt(uMeas.*uMeas+vMeas.*vMeas);
temp=EleSizeIndicator;

Ind=speed>10 & GF.node>0.9 ; temp(Ind)=2*CtrlVar.MeshSizeFastFlow;
Ind=speed>50 & GF.node>0.9 ; temp(Ind)=CtrlVar.MeshSizeFastFlow;

EleSizeIndicator=min(temp,EleSizeIndicator);

%[FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,EleSizeIndicator,CtrlVar);




end
