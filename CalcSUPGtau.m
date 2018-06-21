function [tau1,tau2,ECN,K]=CalcSUPGtau(CtrlVar,MUA,u,v,dt)


%
%  Calculates nodal based tau values to be used in the SUPG method.
%

l=2*sqrt(TriAreaFE(MUA.coordinates,MUA.connectivity));
[M,ElePerNode] = Ele2Nodes(MUA.connectivity,MUA.Nnodes);
l=M*l; 
speed=sqrt(u.*u+v.*v);
ECN=speed.*dt./l;

K=coth(ECN)-1./ECN;  % (1/ECN+ECN/3+..) -1/ECN=ECN/3  if ECN->0 
% turns out the expression for K starts to suffer from numerical errors for ECN < 1e-6
% x=logspace(-10,-5); figure ; semilogx(x,coth(x)-1./x,'r') ; hold on ; semilogx(x, x/3,'b')
I=ECN < 1e-6 ; K(I)=ECN(I)/3 ;  % replaced by the Taylor expansion

tau1=K.*l./speed/2;

% And now I must consider the possibility that speed is zero, in which case
% the above expression fails and must be replaced by the correct limit which is
% tau1 -> dt/6 as speed -> 0
I=speed<100*eps ; tau1(I)=dt/6;


tau2=NaN;

end