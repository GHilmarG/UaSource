

%%
%
% Under development, not for general use
%
%
%
%
%%

% Method options:
%   
%    "Fh + <l , hMeas>"
%    "(h-hmeas) P (h-hmeas) / 2 + <l , Fh>"
%    "(h-hmeas) P (h-hmeas) / 2 + (h-hprior) Q (h-hprior) / 2 + (K h - b) M (K h - b) / 2"
%
%
%

CtrlVar=Ua2D_DefaultParameters(); 

[UserVar,CtrlVar,MUA,F,BCs,Priors,Meas,htrue,kIso,kAlong,kCross,Method]=DefineDataForThicknessInversion(CtrlVar) ; 



PlotDataForThicknessInversion(UserVar,CtrlVar,MUA,F,BCs,Priors,Meas,htrue,kIso,kAlong,kCross,Method) ;


[UserVar,hest,lambda]=hEquation(UserVar,CtrlVar,MUA,F,BCs,kIso,kAlong,kCross,Method,Priors,Meas); 


PlotResultsFromThicknessInversion(CtrlVar,MUA,F,BCs,Priors,Meas,Method,hest,htrue,kIso,kAlong,kCross);






