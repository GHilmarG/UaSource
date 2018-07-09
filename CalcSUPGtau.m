function [tau,tau1,tau2,taus,taut,ECN,K]=CalcSUPGtau(CtrlVar,MUA,u,v,dt)


%
%  Calculates nodal based tau values to be used in the SUPG method.
%

l=sqrt(2*TriAreaFE(MUA.coordinates,MUA.connectivity));
[M,ElePerNode] = Ele2Nodes(MUA.connectivity,MUA.Nnodes);
l=M*l;
speed=sqrt(u.*u+v.*v);

%dt=l./speed;
dt=dt+zeros(MUA.Nnodes,1);

ECN=speed.*dt./l;

K=coth(ECN)-1./ECN;  % (1/ECN+ECN/3+..) -1/ECN=ECN/3  if ECN->0
% turns out the expression for K starts to suffer from numerical errors for ECN < 1e-6
% x=logspace(-10,-5); figure ; semilogx(x,coth(x)-1./x,'r') ; hold on ; semilogx(x, x/3,'b')
I=ECN < 1e-6 ; K(I)=ECN(I)/3 ;  % replaced by the Taylor expansion

tau1=K.*l./speed/2;

% And now I must consider the possibility that speed is zero, in which case
% the above expression fails and must be replaced by the correct limit which is
% tau1 -> dt/6 as speed -> 0
I=speed<100*eps ; tau1(I)=dt(I)/6;


%%
taut=dt/2+eps;
taus=0.5*l./(speed+100*eps);  % Now this must go down to zero gracefully...
tau2=1./(1./taut+1./taus);

switch CtrlVar.Tracer.SUPG.tau
    
    case 'tau1'
        tau=tau1;
    case 'tau2'
        tau=tau2;
    case 'taus'
        tau=taus;
    case 'taut'
        tau=taut;
    otherwise
        error('in CalcSUPGtau case not found')
end

%fprintf('tau1=%f \t tau2=%f \t taus=%f \t taut=%f \n ',mean(tau1),mean(tau2),mean(taus),mean(taut))


%%

end