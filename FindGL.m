function [GLx,GLy,GLxUpper,GLyUpper,GLxLower,GLyLower] = FindGL(DTxy,gfnod,CtrlVar,xBoundary,yBoundary)

    
    % calculates an approximate GL boundary by finding the coordinates of a contourline around all non-floting areas
    % 
        
    F=TriScatteredInterp(DTxy,double(gfnod));
    Resolution=CtrlVar.GLresolutionWhenPlotting;    
    xgrid= min(DTxy.X(:,1)):Resolution:max(DTxy.X(:,1)) ;  
    ygrid= min(DTxy.X(:,2)):Resolution:max(DTxy.X(:,2)) ;  
    
    [X,Y]=meshgrid(xgrid,ygrid);
    
    F.Method='natural';
    grounded=F(X,Y);
    
    % force all locatons outside of meshed domain to be treated as not grounded
    if nargin>3


        
        [In,On] = inpolygon(X,Y,xBoundary,yBoundary);
        
        grounded(In==0)=0;
    end
    
    
    grounded(isnan(grounded))=0;
    
    
    [GLx,GLy]=getcon(X,Y,grounded,0.5);
    if nargout> 3
        [GLxUpper,GLyUpper]=getcon(X,Y,grounded,0.95);
        [GLxLower,GLyLower]=getcon(X,Y,grounded,0.05);
    end
    
    %ind=~isnan(GLx); % don't do this because the NaNs indicate jump to  different grounding lines 
    %GLx=GLx(ind) ; GLy=GLy(ind); % however with the NaNs 2011a crashes!
    %%
end

