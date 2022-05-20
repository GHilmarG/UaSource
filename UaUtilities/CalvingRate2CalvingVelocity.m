function [cx,cy]=CalvingRate2CalvingVelocity(CtrlVar,MUA,F,LSF,c)

% Calculates the calving velocity as a function of the level-set function (LSF) and the calving rate (c)
%
%

[dLSFdx,dLSFdy]=calcFEderivativesMUA(LSF,MUA,CtrlVar);
[nx,ny]=ProjectFintOntoNodes(MUA,dLSFdx,dLSFdy);

Norm=vecnorm([nx ny],2,2);
nx=nx./(Norm+eps) ;
ny=ny./(Norm+eps) ;


cx=c.*nx ; cy=c.*ny;


end

