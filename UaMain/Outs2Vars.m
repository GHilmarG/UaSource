function [CtrlVar,MUA,time,s,b,h,S,B,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,GF,BCs,l]=Outs2Vars(Outs)

%
% [CtrlVar,MUA,time,s,b,h,S,B,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,GF,BCs,l]=Outs2Vars(Outs)
%
%

if isfield(Outs,'s')
    s=Outs.s;
else
    s=[];
end

b=Outs.b;
h=s-b;
S=Outs.S;
B=Outs.B;
ub=Outs.ub;
vb=Outs.vb;
ud=Outs.ud;
vd=Outs.vd;
dhdt=Outs.dhdt;
dsdt=[];
dbdt=[];
C=Outs.C;
AGlen=Outs.AGlen;
m=Outs.m;
n=Outs.n;
rho=Outs.rho;
rhow=Outs.rhow;
g=Outs.g;
as=Outs.as;
ab=Outs.ab;
GF=Outs.GF;
l=Outs.l;
time=Outs.CtrlVar.time;
CtrlVar=Outs.CtrlVar;
MUA=Outs.MUA;
BCs=Outs.BCs;



end


