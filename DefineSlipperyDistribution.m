
function [UserVar,C,m,q,muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

m=3;
C0=3.16e6^(-m)*1000^m*365.2422*24*60*60;

C=C0+zeros(MUA.Nnodes,1);


q=[] ; muk=[] ; % only needed for Coulomb friction type sliding law.

end
