
function ghg_arrow(x,y,vx,vy,velscale,headscale,sharp,head,col,lw,io,xyratio)
                   1 2  3  4     5        6       7     8   9  10 11  12
[nx,ny]=size(x);

if nx ==1
    x=x' ; y=y' ;vx=vx' ;vy=vy';
end

if nargin<=11
    xyratio=1;  % this is the data aspect ratio, should be given by daspect,
    % so sometimes this works:  daspect(daspect) ; t=daspect ; xyratio=t(2)/t(1);

    if nargin <= 10
        io=1;
        if nargin <= 9
            lw=1;
            if nargin <= 8
                col="k";
                if nargin <= 7
                    head =1 ;
                    if nargin <= 6
                        sharp=0.3;
                        if nargin <= 5
                            headscale=0.3;
                            if nargin == 4
                                velscale=1 ;
                            end
                        end
                    end
                end
            end
        end
    end
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


