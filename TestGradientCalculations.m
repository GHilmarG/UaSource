
%%
clear all ; close all
load SSS2dRestart.mat

x=coordinates(:,1); y=coordinates(:,2);
figure ; trisurf(TRIxy,x,y,u,'EdgeColor','none')  ; title('u')

nip=1;
[etaInt,xint,yint,exx,eyy,exy,Eint,e]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);

[dfdx,dfdy,xint,yint]=calcFEderivatives(u,coordinates,connectivity,nip,CtrlVar);

Iint=1;
[Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);

Dx=squeeze(Deriv(:,1,:));

% for each element nEle the x derivative is
unod=reshape(u(connectivity,1),Nele,nod);

nEle=100;

Dx(nEle,:)*unod(nEle,:)'

dudx=Dx(nEle,:)*u(connectivity(nEle,:))





%
%








