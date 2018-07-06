function [taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh,taubxo,taubyo,taubxa,taubya] = BasalDrag(CtrlVar,He,delta,h,B,H,rho,rhow,ub,vb,C,m,uo,vo,Co,mo,ua,va,Ca,ma)


%%
% Returns basal drag and the directional derivatives of basal drag with respect to u,
% v and h.
%
% taubx     :  x component of basal drag
% tauby     :  y component of basal drag
% dtaubxdu  :  derivative of x component of basal drag with respect to u
% dtaubxdv  :  derivative of x component of basal drag with respect to v
% dtaubxdh  :  derivative of x component of basal drag with respect to h
% dtaubydu  :  derivative of y component of basal drag with respect to u
% dtaubydv  :  derivative of y component of basal drag with respect to v
% dtaubydh  :  derivative of y component of basal drag with respect to h
%
%   [ Delta tbx ] =  M    [Delta u]
%   [ Delta tby ]         [Delta v]
%                         [Delta h]
%
% where 
%
%  M = [ dtaubxdu   dtaubxdv   dtaubxdh ]
%      [ dtaubydu   dtaubydv   dtaubydh ]
%
%  Regularisation parameters:
%
%  CtrlVar.Czero
%  CtrlVar.SpeedZero
%  CtrlVar.HeZero
%
%
%
%      He = HeavisideApprox(kH,h-hf,CtrlVar.Hh0);
%   delta = DiracDelta(kH,h-hf,CtrlVar.Hh0);
%      hf = rhow*H./rho;
%       H = S-B
%


%% Basal drag term : ice
% this drag term is zero if the velocities are zero.

beta2i=(C+CtrlVar.Czero).^(-1./m).*(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;

He=He+CtrlVar.HeZero ; % Regularisation 

taubxi=He.*beta2i.*ub; % this is the straightforward (linear) expression for basal stress
taubyi=He.*beta2i.*vb;


% Dbeta2i is zero for m=1. 
Dbeta2i=(1./m-1).*(C+CtrlVar.Czero).^(-1./m).*(ub.^2+vb.^2+CtrlVar.SpeedZero^2).^((1-3*m)./(2*m));

dtaubxdui=He.*(beta2i+Dbeta2i.*ub.*ub);
dtaubydvi=He.*(beta2i+Dbeta2i.*vb.*vb);

dtaubxdvi=He.*(Dbeta2i.*ub.*vb);
dtaubydui=dtaubxdvi;

dtaubxdhi=delta.*beta2i.*ub ;
dtaubydhi=delta.*beta2i.*vb ;



%% Sea ice drag term : ocean

if CtrlVar.IncludeMelangeModelPhysics
    
    
    U=ub-uo;
    V=vb-vo;
    
    
    beta2o=(Co+CtrlVar.Czero).^(-1./mo).*(sqrt(U.*U+V.*V+CtrlVar.SpeedZero^2)).^(1./mo-1) ;
    
    taubxo=(1-He).*beta2o.*U;
    taubyo=(1-He).*beta2o.*V;
    
    Dbeta2o=(1./mo-1).*(Co+CtrlVar.Czero).^(-1./mo).*(U.^2+V.^2+CtrlVar.SpeedZero^2).^((1-3*mo)./(2*mo));
    
    dtaubxduo=(1-He).*(beta2o+Dbeta2o.*U.*U);
    dtaubxdvo=(1-He).*(Dbeta2o.*U.*V);
    
    dtaubyduo=dtaubxdvo;
    dtaubydvo=(1-He).*(beta2o+Dbeta2o.*V.*V);
    
    dtaubxdho=-delta.*beta2o.*U ;
    dtaubydho=-delta.*beta2o.*V ;
    
    
    
    
    U=ub-ua;
    V=vb-va;
    
    
    beta2a=(Ca+CtrlVar.Czero).^(-1./ma).*(sqrt(U.*U+V.*V+CtrlVar.SpeedZero^2)).^(1./ma-1) ;
    
    taubxa=(1-He).*beta2a.*U;
    taubya=(1-He).*beta2a.*V;
    
    Dbeta2a=(1./ma-1).*(Ca+CtrlVar.Czero).^(-1./ma).*(U.^2+V.^2+CtrlVar.SpeedZero^2).^((1-3*ma)./(2*ma));
    
    dtaubxdua=(1-He).*(beta2a+Dbeta2a.*U.*U);
    dtaubxdva=(1-He).*(Dbeta2a.*U.*V);
    
    dtaubydua=dtaubxdva;
    dtaubydva=(1-He).*(beta2a+Dbeta2a.*V.*V);
    
    dtaubxdha=-delta.*beta2a.*U ;
    dtaubydha=-delta.*beta2a.*V ;
    
    
else
    
    taubxo=0;
    taubyo=0;
    
    dtaubxduo=0;
    dtaubxdvo=0;
    
    dtaubyduo=0;
    dtaubydvo=0;
    
    dtaubxdho=0;
    dtaubydho=0;
    
    taubxa=0;
    taubya=0;
    
    dtaubxdua=0;
    dtaubxdva=0;
    
    dtaubydua=0;
    dtaubydva=0;
    
    dtaubxdha=0;
    dtaubydha=0;
end


%% Add up

taubx=taubxi+taubxo+taubxa ;
tauby=taubyi+taubyo+taubya ;

dtaubxdu=dtaubxdui+dtaubxduo+dtaubxdua;
dtaubxdv=dtaubxdvi+dtaubxdvo+dtaubxdva;

dtaubydu=dtaubydui+dtaubyduo+dtaubydua;
dtaubydv=dtaubydvi+dtaubydvo+dtaubydva;

dtaubxdh=dtaubxdhi+dtaubxdho+dtaubxdha;
dtaubydh=dtaubydhi+dtaubydho+dtaubydha;



end