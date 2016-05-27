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
%
%
% Another option of just plotting the mesh boundaries is simply to do: 
% figure ; plot(MUA.coordinates(MUA.Boundary.Edges,1)/CtrlVar.PlotXYscale, MUA.coordinates(MUA.Boundary.Edges,2)/CtrlVar.PlotXYscale, 'k', 'LineWidth',2) ;
% in which case the `PlotBoundary' m-file is not used.
%
%
%%

PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,varargin{:})


end