function [LakeNodes,OceanNodes] = LakeOrOcean3(CtrlVar,MUA,GF)

%%
%
%   [LakeNodes,OceanNodes]=LakeOrOcean3(CtrlVar,MUA,GF)
%
% When calling this to to apply melt, the following syntax is recommended
% to ensure that melt is applied correctly:
%
%   [LakeNodes,OceanNodes]=LakeOrOcean3(CtrlVar,MUA,GF)
%   ab(~OceanNodes) = 0;
%
% This script is designed to be used in conjunction with DefineMassBalance
% to only assign melt to nodes that should be melted (OceanNodes).
% Note that this script does not robustly identify all possible lakes in a 
% domain, since it only considers nodes strictly downstream of the grounding
% line as floating. Thus, floating nodes with an edge that crosses the
% grounding line, which are not considered floating, will also not be
% considered as lakes. In this way, very small isolated patches of floating
% nodes will neither be considered lakes nor ocean.
%
% Also consider using: LakeOrOcean.m , which uses an alternative approach for the problem.
%
%
% Currently, Úa users are split into LakeOrOcean.m and the LakeOrOcean3.m camps.
% The LakeOrOcean3.m approach is to consider lake being a lake if it is enclosed
% by  grounded ice. The LakeOrOcean.m approach is to identify the longest
% grounding line and consider any floating areas upstream of that grounding line
% to be lakes and all other floating areas a part of the ocean. Both of these
% approached can fail. However, arguably the LakeOcean3.m definition of a lake
% is more likely to be generally accepted by members of a typical university
% geography department.
%
%%

GF = IceSheetIceShelves(CtrlVar,MUA,GF);
EleSubset = GF.ElementsDownstreamOfGroundingLines;

FloatingBoundaryNodes=MUA.Boundary.Nodes(GF.node(MUA.Boundary.Nodes)<0.5);

TRI=MUA.connectivity(EleSubset,:) ;

% create undirected graph
G=graph(TRI,TRI(:,[2 3 1]));
% calculate the connected components of the graph
bins=conncomp(G) ;

% initialise arrays
Nnum = zeros(MUA.Nnodes,1);
Nnum(FloatingBoundaryNodes) = 1;
LakeNodes = GF.NodesDownstreamOfGroundingLines;

% loop through ocean boundary nodes until each one has been checked for
% connected floating nodes, once this is done for all boundary nodes the
% only floating nodes left should be lakes
while sum(Nnum)>0
    NodeSeed = find(Nnum,1,'first');
    ID=bins(NodeSeed) ;
    % list of all connected nodes to this ocean boundary node
    nodes=find(bins==ID);
    % remove these from the Lakes list
    LakeNodes(nodes) = 0;
    % also remove these from the list of boundary nodes to save time where
    % one ice shelf has multiple nodes on the ocean boundary
    Nnum(nodes) = 0;
end

OceanNodes = GF.NodesDownstreamOfGroundingLines & ~LakeNodes;
