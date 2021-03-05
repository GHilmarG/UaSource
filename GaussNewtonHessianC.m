function [ddIdCC,UserVar]=GaussNewtonHessianC(UserVar,CtrlVar,MUA,dIdC,F,Meas)


% Gauss Newton estimate of ddIdCC

%
% Straightforward differentaion of I=I(u(C)) = 0.5 ((u-umeas)/uerr)^2
% 
%  dI/dC= dI/du  du/dC   = (u-umes) uerr^{-2}    du/dC
% 
%  ddIdCC=  uerr^{-2}   (du/dC)^2   + second-order derivative with respect to u that we ignore
% 
%
 
uErr=sqrt(spdiags(Meas.usCov));
vErr=sqrt(spdiags(Meas.vsCov));

speed=sqrt(F.ub.*F.ub+F.vb.*F.vb) ;
uMin=max(CtrlVar.SpeedZero,0.01*max(speed)); 

 
gx= uErr./(abs(F.ub-Meas.us)+uMin); 
gy= vErr./(abs(F.vb-Meas.vs)+uMin); 

dIdC=dIdC(:) ; gx=gx(:) ; gy=gy(:) ;


f=(gx.*gx+gy.*gy).*dIdC.*dIdC ;



[UserVar,ddIdCC]=FE_outer_product(UserVar,CtrlVar,MUA,f) ; 

ddIdCC=(ddIdCC+ddIdCC')/2 ; 


end