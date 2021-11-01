function [LakeNodes,OceanNodes,LakeElements,OceanElements] = LakeOrOcean3(CtrlVar,MUA,GF,OceanBoundaryNodes,NodesDownstreamOfGroundingLines)

narginchk(3,5)
nargoutchk(2,4)

%%
%
%   [LakeNodes,OceanNodes] = LakeOrOcean3(CtrlVar,MUA,GF,OceanBoundaryNodes)
%
% OceanBoundaryNodes is an optional input that allows the user to specify
% the node numbers of nodes in their domain that are connected to the
% ocean, if this input is not provided the code assumes that all floating
% nodes on the model boundary are connected to the ocean. This is usually a
% good guess but may break down in certain situations, particularly if
% there are holes in the model mesh.
%
% When calling this to apply melt, the following syntax is recommended
% to ensure that melt is applied correctly:
%
%   [LakeNodes,OceanNodes]=LakeOrOcean3(CtrlVar,MUA,GF)
%   ab(~OceanNodes) = 0;
%
% This script is designed to be used in conjunction with DefineMassBalance
% to only assign melt to nodes that should be melted (OceanNodes). 
%
% Note that this script does not robustly identify all possible lakes in a 
% domain, since it only considers nodes strictly downstream of the grounding
% line as floating. Thus, floating nodes with an edge that crosses the
% grounding line, which are not considered floating, will also not be
% considered as lakes. In this way, very small isolated patches of floating
% nodes will neither be considered lakes nor ocean.
%
% Also consider using: LakeOrOcean.m , which uses an alternative approach for the problem.
% ... but is painfully slow and will fail to correctly identify lakes
% within grounded islands
%
% Currently, Úa users are split into LakeOrOcean.m and the LakeOrOcean3.m camps.
% The author of the LakeOrOcean.m prefers using LakeOrOcean3.m
%
%
%  
%   load PIG-TWG-RestartFile.mat ; CtrlVar=CtrlVarInRestartFile;
%   [LakeNodes,OceanNodes,LakeElements,OceanElements]=LakeOrOcean3(CtrlVar,MUA,F.GF) ;
%
%   FindOrCreateFigure("LakeOrOcean3") ; 
%   CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
%   PlotMuaMesh(CtrlVar,MUA) ;
%   hold on ; plot(MUA.coordinates(OceanNodes,1)/CtrlVar.PlotXYscale,MUA.coordinates(OceanNodes,2)/CtrlVar.PlotXYscale,'ob') ;
%   PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],color='r') ;
%
%%

GF = IceSheetIceShelves(CtrlVar,MUA,GF);

% TestIng
if nargin<5 || isempty(NodesDownstreamOfGroundingLines)
    
    NodesDownstreamOfGroundingLines=GF.NodesDownstreamOfGroundingLines;
else
    if isstring(NodesDownstreamOfGroundingLines)
        if contains(NodesDownstreamOfGroundingLines,"Strickt")
            NodesDownstreamOfGroundingLines=GF.NodesDownstreamOfGroundingLines;
        elseif contains(NodesDownstreamOfGroundingLines,"Relaxed")
            NodesDownstreamOfGroundingLines=GF.node < 0.5 ; 
        end
    end
end

    
% if the user does not provide an OceanBoundaryNodes vector as input,
% assume that floating nodes on the Mesh Boundary are ocean nodes
if nargin<4 || isempty(OceanBoundaryNodes)
    OceanBoundaryNodes=MUA.Boundary.Nodes(GF.node(MUA.Boundary.Nodes)<0.5);
end

% the graph will be comprised only of fully floating elements
EleSubset = GF.ElementsDownstreamOfGroundingLines;

% find the node numbers of all of these fully floating elements
NodeSubset = unique(MUA.connectivity(EleSubset,:));


% don't include floating boundary nodes that are not part of fully floating
% elements as these will not be a part of the graph network
FloatingSubset = intersect(NodeSubset,OceanBoundaryNodes);

TRI=MUA.connectivity(EleSubset,:) ;

% create undirected graph
G=graph(TRI,TRI(:,[2 3 1]));
% calculate the connected components of the graph
bins=conncomp(G) ;

% initialise arrays
Nnum = zeros(MUA.Nnodes,1);

Nnum(FloatingSubset) = 1;
LakeNodes = NodesDownstreamOfGroundingLines;

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

OceanNodes = NodesDownstreamOfGroundingLines & ~LakeNodes;

if nargout > 2
    LakeElements=AllElementsContainingGivenNodes(MUA.connectivity,find(LakeNodes)) ;
end

if nargout > 3
    OceanElements=AllElementsContainingGivenNodes(MUA.connectivity,find(OceanNodes)) ;
end


