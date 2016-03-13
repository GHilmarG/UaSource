function PlotMuaMesh(CtrlVar,MUA,ElementList)

% PlotMuaMesh(CtrlVar,MUA,ElementList)
%
% Examples:
% figure ; PlotMuaMesh([],MUA);
%
% CtrlVar.NodeColor='r';
% figure ; PlotMuaMesh(CtrlVar,MUA,1:100);
%
% CtrlVar.PlotLabels=0;
% CtrlVar.PlotXYscale=1000;
% figure ; PlotMuaMesh(CtrlVar,MUA,1:100);
%

if nargin<3
    
    ElementList=1:MUA.Nele;
end

PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar,ElementList)


end