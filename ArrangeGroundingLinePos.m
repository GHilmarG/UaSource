function [xGL,yGL,GLindex] = ArrangeGroundingLinePos(CtrlVar,GLgeo,NGLmax)
    
    %
    %  Rearanges individual grounding lines found within elements 
    %  to form continous grounding lines. On output xGL and yGL 
    %  are vectors of grounding line (x,y) values with NaN seperating
    %  different grounding lines.
    %
    % First call [GLgeo,GLinfo]=GLgeometry(connectivity,coordinates,GF,CtrlVar)
    % to get GLgeo which is a required input for ArrangeGroundingLinePos
    %
    %  Grounding lines are arranged in order of number of GL points,
    %  with the longest grounding lines first
    %
    %  GLindex(I) gives the index of the first GL point in (xGL,yGL) 
    %  of the Ith grounding line.
    %
    %  numel(GLindex) is the number of returned grounding lines
    %
    % NGLmax is the maximum number of grounding lines returned
    % numel(GLindex) is always equal or smaller than NGLmax
    %
    % To find only the first longest GL and to plot it: 
    % [GLgeo,GLinfo]=GLgeometry(connectivity,coordinates,GF,CtrlVar)
    % [xGL,yGL,GLindex] = ArrangeGroundingLinePos(CtrlVar,GLgeo,1)
    % plot(xGL,yGL)
    %
    %
    %
    
    if isempty(GLgeo) || isempty(~isnan(GLgeo(:,1)))
        xGL=[]; yGL=[] ; GLindex=[];
        return
    end
    
    
    
    if nargin<3
        NGLmax=inf;
    end
    
    
    xa=GLgeo(:,3) ;  xb=GLgeo(:,4) ;
    ya=GLgeo(:,5) ;  yb=GLgeo(:,6) ;
    
    
    [xGL,yGL,GLindex]=LineUpEdges(CtrlVar,xa,xb,ya,yb,NGLmax);

    
end
