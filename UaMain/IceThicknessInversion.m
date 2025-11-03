

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
UserVar=[]; 

[UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar) ; 

[UserVar,CtrlVar,MUA,F,BCs,Priors,Meas,htrue]=DefineDataForThicknessInversion(UserVar,CtrlVar) ; 



PlotDataForThicknessInversion(UserVar,CtrlVar,MUA,F,BCs,Priors,Meas,htrue) ;


[UserVar,hest,lambda]=hEquation(UserVar,CtrlVar,MUA,F,BCs,Priors,Meas); 


PlotResultsFromThicknessInversion(CtrlVar,MUA,F,BCs,Priors,Meas,hest,htrue);






