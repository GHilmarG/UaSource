function B=MismBed(x,y)


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