function [ub,vb,ud,vd,uo,vo,s,b,h,S,B,AGlen,C,m,n,rho,rhow,Co,mo,Ca,ma,as,ab,dasdh,dabdh,dhdt,dsdt,dbdt,dubdt,dvbdt,duddt,dvddt,g,alpha]=UaFields2Vars(F)




ub=F.ub;
vb=F.vb;

ud=F.ud;
vd=F.vd;

uo=F.uo;
vo=F.vo;

s=F.s;
b=F.b;
h=F.h;
S=F.S;
B=F.B;
rho=F.rho;
rhow=F.rhow;

AGlen=F.AGlen;
n=F.n;

C=F.C;
m=F.m;



Co=F.Co;
mo=F.mo;

Ca=F.Ca;
ma=F.ma;

as=F.as;
ab=F.ab;

dasdh=F.dasdh;
dabdh=F.dabdh;

dhdt=F.dhdt;
dsdt=F.dsdt;
dbdt=F.dbdt;

dubdt=F.dubdt;
dvbdt=F.dvbdt;

duddt=F.duddt;
dvddt=F.dvddt;


g=F.g;
alpha=F.alpha;





end
