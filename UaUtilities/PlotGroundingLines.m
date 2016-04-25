function [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,varargin)

% Plots grounding lines
% 
% [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,varargin)
%
% Examples:
%
% PlotGroundingLines(CtrlVar,MUA,GF);
% Plots grounding lines for a given FE-mesh MUA and the floating mask GF.
%
% [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'r')
% Plots grounding lines in red color.
%
% [xGL,yGL,GLgeo]=PlotGroundingLines([],[],[],GLgeo,[],[],'r');
% Plots previously calculated grounding lines.
%
% [xGL,yGL,GLgeo]=PlotGroundingLines([],[],[],[],xGL,yGL,'r'); 
% Plots previously calculated grounding lines.
% 
% Repeated use:
% GLgeo=[] ; xGL=[], yGL=[];
% [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
% [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
%
%
% If xGL and yGL are empty on input, they are calculated from MUA and GLgeo.
% If GLgeo, xGL and yGL are empty on input, then GLgeo is first calculated from from GF and MUA, 
% and then xGL and yGL are calculated from GLgeo.
%
% If GLgeo, xGL and yGL are all empty, then MUA and GF can not be empty.
%
% If xGL and yGL are not empty, then all other fields can be left empty.
%
% varagin is passed on to plot
%
% 
%  Note: If MUA and GF have not changed from a previous call to PlotGroundingLines, always give xGL and yGL from 
%  the previous call as inputs to all following calls as doing so saves time when plotting complicated 
%  grounding lines. In an initial call, where MUA and GF are available, do for example:
%  GLgeo=[] ; xGL=[] ; yGL=[];
% [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL); % here GLgeo, xGL and yGL are all empty
% [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL); % here GLgeo, xGL and yGL are known from the previous call (they will not be recalculated in this call). 
%  
% Once xGL and yGL have been calculated, for example through a previous call to PlotGroundingLines, one can also simply plot the grounding lines as: plot(xGL,yGL) 
%
% To plot individual grounding lines in different colours set CtrlVar.PlotIndividualGLs=1;
%

if isempty(CtrlVar) 
    CtrlVar.XYscale=1; 
    CtrlVar.PlotIndividualGLs=0;
    CtrlVar.PlotGLs=1;
end

if ~isfield(CtrlVar,'PlotXYscale') ; CtrlVar.PlotXYscale=1 ; end
if ~isfield(CtrlVar,'PlotGLs') ; CtrlVar.PlotGLs=1 ; end
if ~isfield(CtrlVar,'PlotIndividualGLs') ; CtrlVar.PlotIndividualGLs=0 ; end


if nargin<4
    GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
    xa=GLgeo(:,3) ;  xb=GLgeo(:,4) ; ya=GLgeo(:,5) ;  yb=GLgeo(:,6) ;
    [xGL,yGL]=LineUpEdges2([],xa,xb,ya,yb);

elseif isempty(xGL) || isempty(yGL)
    
    if isempty(GLgeo)
        
        if isempty(GF)
            error('On input GF can not be empty if xGL, yGL and GLgeo are as well. \n')
        else
            GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
        end
    end
    
    xa=GLgeo(:,3) ;  xb=GLgeo(:,4) ; ya=GLgeo(:,5) ;  yb=GLgeo(:,6) ;
    [xGL,yGL]=LineUpEdges2([],xa,xb,ya,yb);
    
end

if CtrlVar.PlotGLs
    
    plot(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,varargin{:}) ;
    
    
elseif CtrlVar.PlotIndividualGLs
    
    i=0;
    I=find(isnan(xGL)) ;
    I=[1;I(:)];
    col=['b','r','c','g','k','m'];
    %col=['b','r'];
    for ii=1:numel(I)-1
        i=i+1;
        plot(xGL(I(ii):I(ii+1)),yGL(I(ii):I(ii+1)),col(i)) ; axis equal ; hold on ;
        if i==numel(col) ; i=0 ; end
    end
end

ax=gca; ax.DataAspectRatio=[1 1 1];


end