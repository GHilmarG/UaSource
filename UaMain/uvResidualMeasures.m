function [rForce,rWork,D2,Normalisation]=uvResidualMeasures(frhs,grhs,fext0,dub,dvb,dl)

%% Residual measures of the uv system
%
%   [rForce,rWork,D2,Normalisation]=uvResidualMeasures(frhs,grhs,fext0,dub,dvb,dl)
%
% (15 Oct 2026) Moved here, unchanged, from CalcCostFunctionNR.m, so that these measures can also be calculated directly from an
% already assembled residual (eg at the start of a Newton-Raphson iteration, gamma=0), without a further assembly.

Normalisation=fext0'*fext0+1000*eps;
rForce=full([frhs;grhs]'*[frhs;grhs]./Normalisation);

% Newton Decrement
D2=full([frhs;grhs]'*[dub;dvb;dl])  ;
rWork=D2^2 ;

end
