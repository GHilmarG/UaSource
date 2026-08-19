function B=MismBed(x,y)

narginchk(2,2)

% The bedgeometry was defined over a area extending from xmin=0 to xmax=640e3 and ymin=0 to ymax=80e3


%% Here the xy coordinates are scaled, this may, or may not, be what you want...?!
Lx0=640e3; Ly0=80e3; 
xMin=min(x) ; xMax=max(x) ; Lx=xMax-xMin;
yMin=min(y) ; yMax=max(y) ; Ly=yMax-yMin;

x= (x-xMin)/Lx *Lx0;  
y= (y-yMin)/Ly *Ly0;
%%


Ly=80e3;
B0=-150 ;
B2=-728.8 ;
B4=343.91 ;
B6=-50.57;
xbar=300e3;
fc=4e3;
dc=500;
wc=24e3;
Bmax=720;

xtilde=x/xbar;
B=(B0+B2*(xtilde).^2+B4*(xtilde).^4+B6*(xtilde).^6)+...
    dc*(1./(1+exp(-2*(y-Ly/2-wc)/fc))+1./(1+exp(2*(y-Ly/2+wc)/fc)));

B=max(B,-Bmax);

end