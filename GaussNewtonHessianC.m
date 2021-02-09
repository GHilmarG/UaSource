function [ddIdCC,UserVar]=GaussNewtonHessianC(UserVar,CtrlVar,MUA,dIdC,F,Meas)


% Gauss Newton estimate of ddIdCC

%
% Straightforward differentaion of I=I(u(C)) = 0.5 ((u-umeas)/uerr)^2
% 
%  dI/dC= dI/du  du/dC   = (u-umes) uerr^{-2}    du/dC
% 
%  ddIdCC=  uerr^{-2}   (du/dC)^2   + second-order derivative with respect to u that we ignore
% 
% But looking at the relationship between dI/dC and ddI/dCC based on fix point suggests
% H \approx dJ dJ ures^2
%
%

g= (1./uErr)^2 + (1/vErr).^2  ;

if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
   g=log(10)* g ;  
end

% and if inverting for log c I think there is an addionla log(10) factor here as well

f=dIdC(:).*g;




[UserVar,ddIdCC]=FE_outer_product(UserVar,CtrlVar,MUA,f,f) ; 

%Biggest=max(diag(ddIdCC)) ; Smallest=min(diag(ddIdCC)) ; 
%p=0.01*(Biggest-Smallest)+Smallest;
%ddIdCC=ddIdCC+p*MUA.M;


ddIdCC=ddIdCC/MUA.Area; 
ddIdCC=(ddIdCC+ddIdCC')/2 ; 


end