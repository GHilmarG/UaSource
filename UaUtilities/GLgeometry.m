function [GLgeo,GLnodes,GLele]=GLgeometry(connectivity,coordinates,GF,CtrlVar)

%%
% Calculates grounding lines based on the floating/grounding condition at each node.
%
% *Examples:*
%
%       load('MUA-PIG-TWG-Example.mat','MUA','BCs','GF','CtrlVar')
%       [GLgeo,GLnodes,GLele]=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
%
%
% The  grounding line edges can be plotted as:
%
%       figure ; plot(GLgeo(:,[3 4])',GLgeo(:,[5 6])') ; axis equal 
%
% To plot the normals to the grounding line, pointing outwards toward the ocean:
%
%       scale=0.1; H0=1000; hold on ;
%       quiver(GLgeo(:,7)/CtrlVar.PlotXYscale,GLgeo(:,8)/CtrlVar.PlotXYscale,GLgeo(:,9),GLgeo(:,10),scale,'color','r')
%
% CtrlVar is optional, and if given only the CtrlVar.GLsubdivide and
% CtrlVar.GLthreshold fields are needed.
%
%
% CtrlVar.GLsubdivide   :  If CtrlVar.GLsubdivide is set to true, then 6 and 10 nodes elements are subdivided into smaller sub-triangles,
%                          and the grounding line is calculated for each sub-triangle
%
%  CtrlVar.GLthreshold  : A threshold value determining where the GL is located with respect to the values of
%                         flotation mask.  Usually CtrlVar.GLthreshold=0.5
%
%
% Within each element, or sub-triangle in case CtrlVar.GLsubdivide is true, the GL is a straight line
% going from one edge to another edge of the same element/sub-triangle
%
% Output:
% 
%   GLgeo(:,1)                        : list of elements with nodes on both side of the grounding line
%   GLgeo(:,2)                        : edge number
%   GLgeo(:,[3 4])                    : x  coordinates of grounding line position of element GLgeo(:,1) and edge GL(:,2)
%   GLgeo(:,[5 6])                    : y coordinates of grounding line position of element GLgeo(:,1) and edge GL(:,2)
%   GLgeo(:,7)                        : mean x coordinate of grounding line position of element GLgeo(:,1)
%   GLgeo(:,8)                        : mean y coordinate of grounding line position of element GLgeo(:,1)
%   GLgeo(:,9) and GLgeo(:,10)        : (x,y) coordinates of the (outwards) normals to the grounding line
%
%
%   GLnodes : A list of grounded nodes belonging to an element that crosses the grounding line
%      
%
%   GLnodes :  a logical list of nodes belonging to elements crossing grounding lines.
%   GLele :  a logical list of elements crossing grounding lines.
%
% Does a reasonable job, but can only handle one grounding line going through
% each element.
%
% Is fairly well vectorized with running time increasing linearly with problem size,
%
% By default GL position is calculated from GF.node values, using only the corner points and assuming a linear
% variation within element.
%
%
%
% Grounding line of element GLgeo(i,1) goes from edge GL(i,2) to edge rem(GL(i,2),3)+1
%
% If a linked-up list of grounding line positions is required, use:
%       [xGL,yGL,GLindex] = ArrangeGroundingLinePos(CtrlVar,GLgeo,NGLmax)
%
%
% Note on plotting: if the grounding line is plotted as:
%
%            plot(GLgeo(:,[3 4])',GLgeo(:,[5 6])','k')
%
% then matlab creates a new `Line Object' for each and every line segement. For
% grounding lines crossing a large number of elements this can be a very slow
% way of plotting. To speed things up plot the grounding line by first aligning
% up the line segements:
%
%       [xGL,yGL] = ArrangeGroundingLinePos(CtrlVar,GLgeo) ;
%       plot(xGL,yGL,'k')
%
% alternativily, one can also use 
%
%   PlotGroundingLines(CtrlVar,MUA,GF)
%
% to plot the grounding line, which takes care of this issue automatically.
%
% see also: PlotGroundingLines, EleBasedGL
%%

if ~isfield(GF,'node')
    value=GF ;
    clearvars GF
    GF.node=value; 
end

if nargin<4 || isempty(CtrlVar)
    
    CtrlVar.GLthreshold=0.5;
    CtrlVar.GLsubdivide=0;
end


if ~isfield(CtrlVar,'GLthreshold') ; CtrlVar.GLthreshold=0.5; end
if ~isfield(CtrlVar,'GLsubdivide') ; CtrlVar.GLsubdivide=0; end

threshold=CtrlVar.GLthreshold;

if CtrlVar.GLsubdivide
    connectivity=TriFE(connectivity);
end

[Nele,nod]=size(connectivity);
[Nnode,~]=size(coordinates);

switch nod
    case 3
        I=[1 2 3];
    case 6
        I=[1 3 5];
    case 10
        I=[1 4 7];
end


GLele=false(Nele,1);
GLgeo=zeros(Nele,10)+NaN ;  % I'm being a bit wastefull here with memory, but only for a short while as
% this variable is later resized


% going in anti-clockwise direction there can only be one edge for each element for which
% first node has a value less than threshold and the second one a value greater than threshold
%
% I'm assuming all elements have the right orientation, if not the algorithm
% can fail!

for iEdge=1:3
    
    n1=I(rem(iEdge-1,3)+1) ; n2=I(rem(iEdge,3)+1) ; n3=I(rem(iEdge+1,3)+1) ;
    
    GLgeoEdge=(GF.node(connectivity(:,n1)) < threshold) & (GF.node(connectivity(:,n2)) >= threshold ); %
    
    xa=coordinates(connectivity(GLgeoEdge,n1),1) ;  ya=coordinates(connectivity(GLgeoEdge,n1),2) ;
    xb=coordinates(connectivity(GLgeoEdge,n2),1) ;  yb=coordinates(connectivity(GLgeoEdge,n2),2) ;
    xc=coordinates(connectivity(GLgeoEdge,n3),1) ;  yc=coordinates(connectivity(GLgeoEdge,n3),2) ;
    
    gfa=GF.node(connectivity(GLgeoEdge,n1));
    gfb=GF.node(connectivity(GLgeoEdge,n2));
    gfc=GF.node(connectivity(GLgeoEdge,n3));
    
    
    
    GL1x=(threshold-gfa).*(xb-xa)./(gfb-gfa)+xa;
    GL1y=(threshold-gfa).*(yb-ya)./(gfb-gfa)+ya;
    
    
    GL2x=zeros(numel(gfa),1) ;  GL2y=zeros(numel(gfa),1) ;
    
    
    ind1=gfc< threshold ; 
    ind2=gfc>=threshold;
    
    GL2x(ind1)=(threshold-gfb(ind1)).*(xc(ind1)-xb(ind1))./(gfc(ind1)-gfb(ind1))+xb(ind1);  % if gfc<threshold
    GL2x(ind2)=(threshold-gfa(ind2)).*(xc(ind2)-xa(ind2))./(gfc(ind2)-gfa(ind2))+xa(ind2);  % if gfc>threshold
    GL2y(ind1)=(threshold-gfb(ind1)).*(yc(ind1)-yb(ind1))./(gfc(ind1)-gfb(ind1))+yb(ind1);  % if gfc<threshold
    GL2y(ind2)=(threshold-gfa(ind2)).*(yc(ind2)-ya(ind2))./(gfc(ind2)-gfa(ind2))+ya(ind2);  % if gfc>threshold
    
    nx=GL2y-GL1y ; ny=GL1x-GL2x ; temp=sqrt(nx.*nx+ny.*ny+eps); nx=nx./temp ; ny=ny./temp;
    
    
    GLele(GLgeoEdge)=true;
    
    %GLgeo(GLgeoEdge,1)=find(GLgeoEdge) ; 
    GLgeo(GLgeoEdge,2)=iEdge ;
    GLgeo(GLgeoEdge,3)=GL1x ;  
    GLgeo(GLgeoEdge,4)=GL2x ;   
    GLgeo(GLgeoEdge,5)=GL1y ;  
    GLgeo(GLgeoEdge,6)=GL2y ;
    
    GLgeo(GLgeoEdge,7)=(GL1x+GL2x)/2 ;
    GLgeo(GLgeoEdge,8)=(GL1y+GL2y)/2 ;
    
    GLgeo(GLgeoEdge,9)=nx ;
    GLgeo(GLgeoEdge,10)=ny ;
    
    
end

GLgeo(GLele,1)=find(GLele);

% getting rid of all lines in matrix corresponding to elements that do not cross the GL
ind=~isnan(GLgeo(:,1)); 
GLgeo=GLgeo(ind,:);

GLnodes=false(Nnode,1);
GLnodes(connectivity(GLele,:))=true;

end