function F=Vars2UaFields(ub,vb,ud,vd,uo,vo,s,b,h,S,B,AGlen,C,m,n,rho,rhow,Co,mo,Ca,ma,as,ab,dasdh,dabdh,dhdt,dsdt,dbdt,dubdt,dvbdt,duddt,dvddt,g,alpha)


F=UaFields;

F.ub=ub;
F.vb=vb;

F.ud=ud;
F.vd=vd;

F.uo=uo;
F.vo=vo;

F.s=s;
F.b=b;
F.h=h;
F.S=S;
F.B=B;
F.rho=rho;
F.rhow=rhow;

F.AGlen=AGlen;
F.n=n;

F.C=C;
F.m=m;



F.Co=Co;
F.mo=mo;

F.Ca=Ca;
F.ma=ma;

F.as=as;
F.ab=ab;

F.dasdh=dasdh;
F.dabdh=dabdh;

F.dhdt=dhdt;
F.dsdt=dsdt;
F.dbdt=dbdt;

F.dubdt=dubdt;
F.dvbdt=dvbdt;

F.duddt=duddt;
F.dvddt=dvddt;


F.g=g;
F.alpha=alpha;





end
