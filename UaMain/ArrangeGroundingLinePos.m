function [xGL,yGL] = ArrangeGroundingLinePos(CtrlVar,GLgeo,NGLmax)
    
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
    %
    %
    % NGLmax is the maximum number of grounding lines returned
    %
    % To find only the first longest GL and to plot it: 
    % [GLgeo,GLinfo]=GLgeometry(connectivity,coordinates,GF,CtrlVar)
    % [xGL,yGL] = ArrangeGroundingLinePos(CtrlVar,GLgeo,1)
    % plot(xGL,yGL)
    %
    %
    %
    
    if isempty(GLgeo) || isempty(~isnan(GLgeo(:,1)))
        xGL=[]; yGL=[] ; 
        return
    end
    
    
    
    if nargin<3
        NGLmax=inf;
    end
    
    
    xa=GLgeo(:,3) ;  xb=GLgeo(:,4) ;
    ya=GLgeo(:,5) ;  yb=GLgeo(:,6) ;
    
    
    [xGL,yGL]=LineUpEdges2(CtrlVar,xa,xb,ya,yb,NGLmax);

    
end
