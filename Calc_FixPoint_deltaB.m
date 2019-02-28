function dIdB=Calc_FixPoint_deltaB(CtrlVar,MUA,F,Meas)


speed=sqrt(F.ub.*F.ub+F.vb.*F.vb) ;
speedMeas=sqrt(Meas.us.*Meas.us+Meas.vs.*Meas.vs)  ;
[dsdx,dsdy]=calcFEderivativesMUA(F.s,MUA,CtrlVar) ;
slope=sqrt(dsdx.*dsdx+dsdy.*dsdy) ;  % this should be based on measurements
slope=ProjectFintOntoNodes(MUA,slope) ;
tau= F.rho.*F.g.*F.h.*slope ;
dIdB= -(speed-speedMeas).* F.m .* F.C .* tau.^(F.m)./F.h ; % if speed > speedMeas, decreas B

dIdB=dIdB.*F.GF.node ;

% guard agains very large elements in gradient
% I=isoutlier(dIdB,'quartiles','ThresholdFactor',2) ;  dIdB(I)=0 ; 


end