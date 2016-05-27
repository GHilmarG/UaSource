function PlotMuaMesh(CtrlVar,MUA,ElementList)

%%
% PlotMuaMesh(CtrlVar,MUA,ElementList)
%
% Examples:
%
% figure ; PlotMuaMesh(CtrlVar,MUA,ElementList)
% 
% figure ; PlotMuaMesh([],MUA);  % CtrlVar is an optional input
%
% CtrlVar.NodeColor='r';
% figure ; PlotMuaMesh(CtrlVar,MUA,1:100);  % Plot nodes in red
%
% CtrlVar.PlotLabels=0;
% CtrlVar.PlotXYscale=1000;
% figure ; PlotMuaMesh(CtrlVar,MUA,1:100);  % Show only elemetns 1 to 100
%%

if nargin<3
    
    ElementList=1:MUA.Nele;
end

PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar,ElementList)


end