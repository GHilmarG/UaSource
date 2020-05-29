function [LakeNodes,OceanNodes] = LakeOrOcean2(CtrlVar,MUA,GF,GLgeo)
%
% [NodesOcean,NodesLakes]=LakeOrOcean(CtrlVar,MUA,GF, GLgeo)
%
% What is a lake? I would argue it is any group of floating nodes that is
% entirely surrounded by land. 
% The idea in this script is to go through every grounding line polygon 
% that is a closed loop and check whether the nodes inside are floating or
% grounded. 

% This script could conceivably miss lakes if floating nodes appear right 
% on the boundary of the domain since they will not form lakes that are 
% fully enclosed within a grounding line polygon...

if nargin<4 || isempty(GLgeo)
    [GLgeo,~,~]=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
end


x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
[xGL,yGL] = ArrangeGroundingLinePos(CtrlVar,GLgeo);

% I=find(GF.NodesDownstreamOfGroundingLines);
I = find(GF.node<CtrlVar.GLthreshold);

% create index of where this vector starts a new GL
if any(isnan(xGL))
    GL_ind = [0; find(isnan(xGL))];
else
    GL_ind = [0; numel(xGL)];
end

LakeNodes=false(MUA.Nnodes,1);

for ii=1:numel(GL_ind)-1 % loop through every GL polyline
    if xGL(GL_ind(ii)+1)-xGL(GL_ind(ii+1)-1) == 0
        if yGL(GL_ind(ii)+1)-yGL(GL_ind(ii+1)-1) == 0 
            % start and end point are the same... ie. this is an enclosed GL polygon... 
            % so is it a lake or island?
            % if the GL polygon encloses floating nodes, it must be a Lake
            IN = inpoly([x(I) y(I)],[xGL(GL_ind(ii)+1:GL_ind(ii+1)-1) yGL(GL_ind(ii)+1:GL_ind(ii+1)-1)],[],1);
            if sum(IN)>0
                LakeNodes(I(IN)) = true;
            end
        end
    end
end
    
OceanNodes = false(MUA.Nnodes,1);
OceanNodes(I) = true;
OceanNodes(LakeNodes) = false;


if CtrlVar.doplots && CtrlVar.PlotOceanLakeNodes
    
    figure
    hold off
    
    plot(x(OceanNodes)/CtrlVar.PlotXYscale,y(OceanNodes)/CtrlVar.PlotXYscale,'og','MarkerFaceColor','g','DisplayName','Ocean Nodes') ; hold on
    plot(x(LakeNodes)/CtrlVar.PlotXYscale,y(LakeNodes)/CtrlVar.PlotXYscale,'or','MarkerFaceColor','r','DisplayName','Lake Nodes')
    legend

    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,[],[],'r');

    axis equal tight
    hold on
    title('Ocean/Lake nodes')
    
    
end