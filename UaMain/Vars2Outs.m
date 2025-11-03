function Outs=Vars2Outs(CtrlVar,MUA,s,b,S,B,h,ub,vb,ud,vd,uo,vo,dhdt,C,AGlen,m,n,rho,rhow,g,as,ab,dasdh,dabdh,GF,BCs,l)

Outs.time=CtrlVar.time;
Outs.dt=CtrlVar.dt;
Outs.s=s;
Outs.b=b;
Outs.h=h;
Outs.S=S;
Outs.B=B;
Outs.ub=ub;
Outs.vb=vb;
Outs.ud=ud;
Outs.vd=vd;
Outs.uo=uo;
Outs.vo=vo;
Outs.dhdt=dhdt;
Outs.C=C;
Outs.AGlen=AGlen;
Outs.m=m;
Outs.n=n;
Outs.rho=rho;
Outs.rhow=rhow;
Outs.g=g;
Outs.as=as;
Outs.ab=ab;
Outs.dasdh=dasdh;
Outs.dabdh=dabdh;
Outs.GF=GF;
Outs.l=l;

Outs.CtrlVar=CtrlVar;
Outs.MUA=MUA;
Outs.BCs=BCs;


end