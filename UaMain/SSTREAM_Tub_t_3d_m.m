function [trans]=SSTREAM_Tub_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m)


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
gamm=tau^(1-m) / (C*m); % can't use gamma because it is a function

ps=m.*gamm.^2+H.*(l.^2.*(4+m)+k.^2.*(1+4.*m)).*gamm.*eta+4.*H.^2.*(k.^2+l.^2).^2.*eta.^2; % better not use psi because it is a standard matlab function

%kappa=l.^2.*(4+m)+k.^2.*(1+4*m);


 p=SSTREAM_Pole_t_3d_m(kx,ky,alpha,H,eta,C,rho,g,m);

if isnan(t)
  expp=0;
  %disp(' steady state option used ')
else
  expp=exp(p*t);
end

%t1=-tau.*(ca.*H.*(k.^2.*U.*(m.*gamm+H.*j2.*eta)-l.^2.*tau)+expp.*(m.*(p-ca.*H.*k.^2.*U).*gamm+H.*(k.^2.*(p-ca.*H.*j2.*U).*eta+l.^2.*(4.*p.*eta+ca.*tau))));
%t2=H*p.*(m.*gamm^2+H.*eta.*(4.*H.*j2.^2.*eta+gamm.*kappa));

t1=-tau.*(ca.*H.*(k.^2.*U.*(m.*gamm+H.*j2.*eta)-l.^2.*tau)+expp.*(m.*(p-ca.*H.*k.^2.*U).*gamm+H.*(k.^2.*(p-ca.*H.*j2.*U).*eta+l.^2.*(4.*p.*eta+ca.*tau))));

t2=H.*p.*ps;

trans=t1./t2;

return

end



