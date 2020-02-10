function [UserVarInRestartFile,CtrlVarInRestartFile,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo]=PlotInversion(InverseRestartFile)


%%
%  
%   PlotInversion(InverseRestartFile)
%
% Plots results in an inverse restart file.
%
%%
if nargin==0
    
    % I'm assuming here that the name of the inverse restart file starts with the letters IR
    InverseRestartFile=uigetfile(...
        "IR*.mat",...
        "Select inverse restart file to plot");
    
end

load(InverseRestartFile,'UserVarInRestartFile','CtrlVarInRestartFile','MUA','BCs','F','l','InvStartValues','InvFinalValues','Priors','Meas','BCsAdjoint','RunInfo');

PlotResultsFromInversion(UserVarInRestartFile,CtrlVarInRestartFile,MUA,BCs,F,l,F.GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);

end