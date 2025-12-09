function [trans]=SSTREAM_Tvb_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m)


% time-dependent ratio between forward surface velocity component and
% basal slipperiness Fourier space in dimensional units

% for non-dimensional a la gudmundsson 2003
% put eta=1/2, H0=1,rho*g=1/sin(alpha), lambda in units of mean ice thickness
% Cnondimensional=2 eta/H0 Cdimesional, and then
% U0nondimensional=Cnondimensional 

[k,l] = kxky2kl(kx,ky);

j2=k.^2+l.^2;
ca=cot(alpha);
tau=rho*g*sin(alpha)*H;
U=C*tau^m;

%kappa=l.^2.*(4+m)+k.^2.*(1+4*m);
%nu=gamm+H*j2*eta;

gamm=tau^(1-m) / (C*m); % can't use gamma because it is a function
ps=m.*gamm.^2+H.*(l.^2.*(4+m)+k.^2.*(1+4.*m)).*gamm.*eta+4.*H.^2.*j2.^2.*eta.^2;

p=SSTREAM_Pole_t_3d_m(kx,ky,alpha,H,eta,C,rho,g,m);



if isnan(t)
  expp=0;
  %disp(' steady state option used ')
else
  expp=exp(p*t);
end

%t1=k.*l.*tau.*(3*expp.*p.*eta+ ca*(expp-1).*(U*nu+tau));
%t2=p.*(m*gamm^2+H*eta*(4*H*j2.^2*eta+gamm*kappa));

t1=k.*l.*tau.*(3.*expp.*p.*eta+ca.*(expp-1).*(U.*(gamm+H.*j2.*eta)+tau));
t2=p.*ps;

trans=t1./t2;

return

end



