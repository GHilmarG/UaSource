
function [OceanNodes,LakeNodes]=LakeOrOcean(CtrlVar,GF,Boundary,connectivity,coordinates)
%
% [NodesOcean,NodesLakes]=LakeOrOcean(CtrlVar,GF,Boundary,connectivity,coordinates)
%
% Tries to determine which floating nodes are part of the ocean and which belong to subglacial lakes
%
%    uses GLgeometry to determine the longest grounding line and assumes that this grounding line
%    represents the ocean boundary. Then determines if nodal point at floating level are inside or
%    outside of this boundary.
%
% May not always work, and anyhow I am not sure if one can always objectivily decide what is an ocean
% and what a lake.
%
% Note: GLgeo calculated is slighly different from the usual way of doing this
%       because here the grounding line needs to be closed.
%%
OceanNodes=[];
LakeNodes=[];


if ~isfield(CtrlVar,'GLthreshold')
    CtrlVar.GLthreshold=0.5;
end

I=find(GF.node<CtrlVar.GLthreshold);   % all floating nodes

if ~isempty(I)
    
    
    GFtemp=GF; GFtemp.node(Boundary.Nodes)=0;  % I need to `close' the grounding line, so I set all boundary nodes to floating status
    GLgeo=GLgeometry(connectivity,coordinates,GFtemp,CtrlVar);
    [xGL,yGL] = ArrangeGroundingLinePos(CtrlVar,GLgeo,1);
    
    x=coordinates(:,1); y=coordinates(:,2);
    %  IN = inpolygon(x(I),y(I),xGL,yGL);  % for some reason this standard matlab routine is much slower than inpoly
    [IN,ON] = inpoly([x(I) y(I)],[xGL yGL],[],1);
    
    % There is a bit of a question here what to do with nodes that are
    % directly on the grounding line. I've here decided to consider them part of
    % the ocean.
    Ind=~IN  | ON ; % ocean nodes are not within grounding line, but can be on it
    
    OceanNodes=I(Ind);
    LakeNodes=I(~Ind);
    
    
    if CtrlVar.doplots && CtrlVar.PlotOceanLakeNodes
        
        figure
        hold off
        
        plot(x(OceanNodes)/CtrlVar.PlotXYscale,y(OceanNodes)/CtrlVar.PlotXYscale,'+g') ; hold on
        plot(x(LakeNodes)/CtrlVar.PlotXYscale,y(LakeNodes)/CtrlVar.PlotXYscale,'or')
        
        
        plot(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,'k','LineWidth',2);
        plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',1);
        %plot(x(Boundary.EdgeCornerNodes)/CtrlVar.PlotXYscale,y(Boundary.EdgeCornerNodes)/CtrlVar.PlotXYscale,'k.-')
        %plot(x(Boundary.Nodes)/CtrlVar.PlotXYscale,y(Boundary.Nodes)/CtrlVar.PlotXYscale,'ro')
        axis equal tight
        hold on
        title('Ocean/Lake nodes')
        legend('Ocean','Lake','Coast','GL')
        
    end
end
end