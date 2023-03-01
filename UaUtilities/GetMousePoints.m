function  [xp,yp]=GetMousePoints(xp,yp,PlotPoints)

% get (x,y) coordinates from mouse clicks
% only selects left-clicks
% returns if right or middle mouse buttons pressed
% if xp yp given as inputs, then appends

N=1000 ;

if nargin <3 
    PlotPoints=1;
end

if nargin>=2 && ~isempty(xp) 
    iCount=numel(xp); 
    xp=[xp(:);zeros(N,1)+NaN] ; yp=[yp(:);zeros(N,1)+NaN];
else 
    xp=zeros(N,1)+NaN ; yp=zeros(N,1)+NaN;
    iCount=1;
end

hold on
button=1;
while button==1
    [xx,yy,button]=ginput(1);
    if button==1
        iCount=iCount+1;
        xp(iCount)=xx ; yp(iCount)=yy;
        if PlotPoints
            plot(xp(iCount),yp(iCount),'+r')
        end
    end
end

I=~isnan(xp) ; xp=xp(I) ; yp=yp(I);


end