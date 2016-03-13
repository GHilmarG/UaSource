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


PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,varargin{:})


end