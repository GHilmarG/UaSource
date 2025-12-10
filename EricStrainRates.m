function [e,exx,eyy,exy,xint,yint,uMeas,vMeas]=EricStrainRates(CtrlVar,MUA,AGlen,n,uMeas,vMeas)


if nargin<4
    AGlen=[];
    n=[];
end

if nargin<5
    
    [uMeas,vMeas]=EricVelocities(CtrlVar,MUA.coordinates);
end

[etaInt,xint,yint,exx,eyy,exy,Eint,e]=calcStrainRatesEtaInt(CtrlVar,MUA,uMeas,vMeas);

%[exxNod,eyyNod,exyNod,eNod]=ProjectFintOntoNodes(MUA,exx,eyy,exy,e)


end