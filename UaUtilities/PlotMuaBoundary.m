function PlotMuaBoundary(CtrlVar,MUA,varargin)

%%
% Plots the boundaries of the MUA mesh.
% PlotMuaBoundary(CtrlVar,MUA,varargin)
%
% varargin is passed over to plot
%
% Examples: 
% PlotMuaBoundary(CtrlVar,MUA)
% PlotMuaBoundary(CtrlVar,MUA,'b')  % plots the boundaries in blue
% PlotMuaBoundary([],MUA,'r','LineWidth',2)

PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,varargin{:})


end