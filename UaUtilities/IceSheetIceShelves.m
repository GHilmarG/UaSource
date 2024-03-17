function [GF,GLgeo,GLnodes,GLele]=IceSheetIceShelves(CtrlVar,MUA,GF,GLgeo,GLnodes,GLele)

%%
% Adds further fields to GF listing elements and nodes upstream,
% downstream, and over grounding lines.
%
%   [GF,GLgeo,GLnodes,GLele]=IceSheetIceShelves(CtrlVar,MUA,GF,GLgeo,GLnodes,GLele)
%
% GLgeo, GLnodes and GLele are optional, but giving them as an input, for example
% from a previous call, speeds things up. (But only do so if GF and MUA has not
% changed between calls!)
%
% All calculated fields are logical variables, ie not nodal numbers. 
%
%
% Element is defined as crossing a grounding line if not all of its nodes are either
% afloat or grounded, as based on the value of the  GL.node floating mask with respect to the
% CtrlVar.GLthreshold value (typically set to 0.5).
%
% Nodes are defined as crossing a grounding line if they are nodes of elements
% that cross a grounding line.
% Example:
%
%   load('MUA-PIG-TWG-Example.mat','MUA','GF','CtrlVar')
%   [GF,GLgeo,GLnodes,GLele]=IceSheetIceShelves(CtrlVar,MUA,GF);
%    
%   figure
%   PlotMuaMesh(CtrlVar,MUA,GF.ElementsUpstreamOfGroundingLines,'color','k')
%   hold on
%   PlotMuaMesh(CtrlVar,MUA,GF.ElementsDownstreamOfGroundingLines,'color','b')
%   PlotMuaMesh(CtrlVar,MUA,GF.ElementsCrossingGroundingLines,'color','r')
% 
%   x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2) ;
%   figure
% 
%   hold on
%   plot(x(GF.NodesDownstreamOfGroundingLines)/CtrlVar.PlotXYscale,y(GF.NodesDownstreamOfGroundingLines)/CtrlVar.PlotXYscale,'.b')
%   hold on
%   plot(x(GF.NodesUpstreamOfGroundingLines)/CtrlVar.PlotXYscale,y(GF.NodesUpstreamOfGroundingLines)/CtrlVar.PlotXYscale,'.k')
%   plot(x(GF.NodesCrossingGroundingLines)/CtrlVar.PlotXYscale,y(GF.NodesCrossingGroundingLines)/CtrlVar.PlotXYscale,'.r')
%   axis equal
%%
%
% SEE ALSO GLgeometry
%
%%

narginchk(3,6)
nargoutchk(1,4)

GF.ele=Nodes2EleMean(MUA.connectivity,GF.node);

if nargin<6 || isempty(GLgeo) || isempty(GLnodes) || isempty(GLele)
    [GLgeo,GLnodes,GLele]=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
end


GF.ElementsCrossingGroundingLines=GLele;
GF.ElementsUpstreamOfGroundingLines=   (GF.ele>CtrlVar.GLthreshold)  & ~GLele;
GF.ElementsDownstreamOfGroundingLines= (GF.ele<CtrlVar.GLthreshold)  & ~GLele;

GF.NodesCrossingGroundingLines=GLnodes;
GF.NodesUpstreamOfGroundingLines=false(MUA.Nnodes,1);
GF.NodesDownstreamOfGroundingLines=false(MUA.Nnodes,1);

GF.NodesUpstreamOfGroundingLines(MUA.connectivity(GF.ElementsUpstreamOfGroundingLines,:))=true;
GF.NodesDownstreamOfGroundingLines(MUA.connectivity(GF.ElementsDownstreamOfGroundingLines,:))=true;

% Now I must make sure that the nodes along the edges of elements crossing the
% grounding line are not included in the up and down stream nodes.
%
GF.NodesDownstreamOfGroundingLines=GF.NodesDownstreamOfGroundingLines & ~GF.NodesCrossingGroundingLines;
GF.NodesUpstreamOfGroundingLines=GF.NodesUpstreamOfGroundingLines & ~GF.NodesCrossingGroundingLines;


return



end