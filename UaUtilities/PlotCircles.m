function PlotCircles(x,y,r,varargin)

    x=x(:) ; y=y(:) ; r=r(:);    
    t=linspace(0,2*pi,25);
    
    for I=1:numel(x)
        
        xx=x(I)+r(I)*cos(t);
        yy=y(I)+r(I)*sin(t);
        plot(xx(:),yy(:),varargin{:})
    end
    
    
end

