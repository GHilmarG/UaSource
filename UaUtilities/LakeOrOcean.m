
function [OceanNodes,LakeNodes,GLgeo,GLnodes,GLele,OceanElements]=LakeOrOcean(CtrlVar,MUA,GF,GLgeo,GLnodes,GLele)

%%
%
%   [OceanNodes,LakeNodes,GLgeo,GLnodes,GLele]=LakeOrOcean(CtrlVar,MUA,GF,GLgeo,GLnodes,GLele)
%
%
%  Only CtrlVar and MUA are required input variables:
%
% Example:
%
%   load PIG-TWG-RestartFile.mat ; CtrlVar=CtrlVarInRestartFile;
%   CtrlVar.doplots=1 ;  CtrlVar.PlotOceanLakeNodes=1;
%   [OceanNodes,LakeNodes,GLgeo,GLnodes,GLele]=LakeOrOcean(CtrlVar,MUA,F.GF) ;
%
%
%
% Tries to determine which floating nodes are part of the ocean and which belong
% to subglacial lakes.
%
% Uses GLgeometry to determine the longest grounding line and assumes that this
% grounding line represents the ocean boundary. Then determines if nodal point
% at floating level are inside or outside of this boundary.
%
% May not always work, and anyhow I am not sure if one can always objectively
% decide what is an ocean and what a lake.
%
% This approach will fail in some circumstances. If, for example, parts of an
% ice-rise go afloat and the resulting subglacial lake does not make a
% connection to the surronding ocean, then this lake will not be identified as a
% lake but classified as a part of the ocean...
%
% Just use this m-file as a starting point towards defining lake/ocean nodes.
%
% Note: GLgeo calculated is slighly different from the usual way of doing this
%       because here the grounding line needs to be closed.
%
% Returns a logical indexing (this was changed from indexing vectors on 20 Dec,
% 2018)
%
%
% Note: An alternative starting point might be:
%
%   IEle=EleFlooding(CtrlVar,MUA,NodeSeed,EleSubset);
%
% Also consider using: LakeOrOcean3.m , which uses an alternative approach for
% the problem.
%
% Currently, Úa users are split into LakeOrOcean.m and the LakeOrOcean3.m camps.
% The LakeOrOcean3.m approach is to consider lake being a lake if it is enclosed
% by  grounded ice. The LakeOrOcean.m approach is to identify the longest
% grounding line and consider any floating areas upstream of that grounding line
% to be lakes and all other floating areas a part of the ocean. Both of these
% approached can fail. ]


% However, arguably the LakeOcean3.m definition of a lake
% is more likely to be generally accepted by members of a typical university
% geography department.
%
%%
OceanNodes=[];
LakeNodes=[];
OceanElements=[]; 


if nargin<4 || isempty(GLgeo) || isempty(GLnodes) || isempty(GLele)
    [GLgeo,GLnodes,GLele]=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
end

if ~isfield(GF,'NodesDownstreamOfGroundingLines')
    [GF,GLgeo,GLnodes,GLele]=IceSheetIceShelves(CtrlVar,MUA,GF,GLgeo,GLnodes,GLele);
end

I=find(GF.NodesDownstreamOfGroundingLines);


if ~isempty(I)
    
    
    GFtemp=GF; GFtemp.node(MUA.Boundary.Nodes)=0;  % I need to `close' the grounding line, so I set all boundary nodes to floating status
    GLgeoMod=GLgeometry(MUA.connectivity,MUA.coordinates,GFtemp,CtrlVar);
    [xGL1,yGL1] = ArrangeGroundingLinePos(CtrlVar,GLgeoMod,1);
    
    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    %  IN = inpolygon(x(I),y(I),xGL,yGL);  % for some reason this standard matlab routine is much slower than inpoly
    [IN,ON] = inpoly2([x(I) y(I)],[xGL1 yGL1],[],1);
    
    % There is a bit of a question here what to do with nodes that are
    % directly on the grounding line. I've here decided to consider them part of
    % the ocean.
    Ind=~IN  | ON ; % ocean nodes are not within grounding line, but can be on it
    
    
    OceanNodes=I(Ind);
    LakeNodes=I(~Ind);
    
    II=false(MUA.Nnodes,1);
    II(OceanNodes)=true;
    OceanNodes=II;
    
    II=false(MUA.Nnodes,1);
    II(LakeNodes)=true;
    LakeNodes=II;
    
    OceanElements=AllElementsContainingGivenNodes(MUA.connectivity,find(OceanNodes)) ; 
    

    if CtrlVar.doplots && CtrlVar.PlotOceanLakeNodes
        
        figure
        hold off
        
        
        P1=plot(x(OceanNodes)/CtrlVar.PlotXYscale,y(OceanNodes)/CtrlVar.PlotXYscale,'og','DisplayName','Ocean Nodes') ; hold on
        P2=plot(x(LakeNodes)/CtrlVar.PlotXYscale,y(LakeNodes)/CtrlVar.PlotXYscale,'or','DisplayName','Lake Nodes');
        
        
        PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar) ; hold on 
        plot(xGL1/CtrlVar.PlotXYscale,yGL1/CtrlVar.PlotXYscale,'k','LineWidth',2);
        [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,[],[],'r');
        %plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',1);
        %plot(x(Boundary.EdgeCornerNodes)/CtrlVar.PlotXYscale,y(Boundary.EdgeCornerNodes)/CtrlVar.PlotXYscale,'k.-')
        %plot(x(Boundary.Nodes)/CtrlVar.PlotXYscale,y(Boundary.Nodes)/CtrlVar.PlotXYscale,'ro')
        axis equal tight
        hold on
        title('Ocean/Lake nodes')
        lg=legend([P1 P2]);
        
    end
end
end