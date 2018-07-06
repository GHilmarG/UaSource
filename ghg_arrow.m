
function ghg_arrow(x,y,vx,vy,velscale,headscale,sharp,head,col,lw,io,xyratio)
 
[nx,ny]=size(x);

if nx ==1
    x=x' ; y=y' ;vx=vx' ;vy=vy';
end
    
    
if nargin == 4
    velscale=1 ;
    headscale=0.3;
    sharp=0.3;
    head=1;
    col='k';
end

if nargin == 5
    headscale=0.3;
    sharp=0.3;
    head=1;
    col='k';
end

if nargin <= 9
    lw=1;
end

if nargin <= 10
    io=1;
end
    
if nargin<=11
    xyratio=1;  % this is the data aspect ratio, should be given by daspect,
                % so sometimes this works:  daspect(daspect) ; t=daspect ; xyratio=t(2)/t(1);
end
vx=vx*velscale/xyratio ; vy=vy*velscale ;



tip_x=x+vx; 
tip_y=y+vy;

X=[x' ; tip_x'];
Y=[y' ; tip_y'];


if io==-1
	tip_x=x;
	tip_y=y;
	
	X=[x'+vx' ; tip_x'];
	Y=[y'+vy' ; tip_y'];
	head=-1;
end

line(X,Y,'color',col,'LineWidth',lw)


ol_x=tip_x-headscale*(head*vx*xyratio+sharp*vy)/xyratio ; ol_y=tip_y+headscale*(-head*vy+sharp*vx*xyratio) ;
X=[tip_x' ; ol_x'] ; Y=[tip_y' ; ol_y'] ;
line(X,Y,'color',col,'LineWidth',lw)

or_x=tip_x+headscale*(-head*vx*xyratio+sharp*vy)/xyratio ; or_y=tip_y+headscale*(-head*vy-sharp*vx*xyratio) ;
X=[tip_x' ; or_x'] ; Y=[tip_y' ; or_y'] ;
line(X,Y,'color',col,'LineWidth',lw)

return


