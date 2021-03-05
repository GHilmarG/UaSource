function [ddIdAA,UserVar]=GaussNewtonHessianA(UserVar,CtrlVar,MUA,dIdA,F,Meas)


% Gauss Newton estimate of ddIdCC

% uErr=sqrt(spdiags(Meas.usCov)); vErr=sqrt(spdiags(Meas.vsCov));

% g=((uErr./F.ub).^2+(vErr./F.vb).^2);

g=1;


f=dIdA(:).*g;


% if dIdA=0 for some i, the H_{ii} will be zero, this will be OK in principle since both H and dIdA are now zero
% 
% ddIdAA=(f*f').*MUA.M ; % f'f will be large!

[UserVar,ddIdAA]=FE_outer_product(UserVar,CtrlVar,MUA,f,f) ; 

Biggest=max(diag(ddIdAA)) ; Smallest=min(diag(ddIdAA)) ; 
p=0.01*(Biggest-Smallest)+Smallest;
ddIdAA=ddIdAA+p*MUA.M;


ddIdAA=ddIdAA/MUA.Area; 
ddIdAA=(ddIdAA+ddIdAA')/2 ; 


end