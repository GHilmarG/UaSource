
function [trans]=SSTREAM_Tsb_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m)


% time-dependent ratio between surface topography and
% bedrock in Fourier space in dimensional units

% for non-dimensional a la gudmundsson 2003
% put eta=1/2, H0=1,rho*g=1/sin(alpha), lambda in units of mean ice thickness
% Cnondimensional=2 eta/H0 Cdimesional, and then
% U0nondimensional=Cnondimensional 


[k,l] = kxky2kl(kx,ky);


j2=k.^2+l.^2;
tau=rho*g*sin(alpha)*H;
U=C*tau^m;
gamm=tau^(1-m) / (C*m); % better not use gamma because it is a standard matlab function
ps=m.*gamm.^2+H.*(l.^2.*(4+m)+k.^2.*(1+4.*m)).*gamm.*eta+4.*H.^2.*(k.^2+l.^2).^2.*eta.^2; % better not use psi because it is a standard matlab function

p=SSTREAM_Pole_t_3d_m(kx,ky,alpha,H,eta,C,rho,g,m);

if isnan(t)
  expp=0;
  %disp(' steady state option used ')
else
  expp=exp(p*t);
end




%t1=(-1i*(-1 + expp).*k.*(U*(gamm + 4*H.*j2*eta) + tau)); %pre Oct 2009 version, same as TC paper
%t2=p.*(gamm + 4*H*j2*eta);




% t1=-1i*(-1 + expp).*k.*(m*gamm*(U*(gamm+H*(4*k.^2+l.^2)*eta)+tau)+H*eta*(k.^2+4*l.^2)*U*gamm+4*H*j2.^2*U*eta+j2*tau);
% t2=p.*(m*gamm^2+H*eta*(4*H*j2.^2*eta+gamm*kappa));


% 07 April 2010
t1=-1i*(expp-1).*k.*(m.*gamm.*(U.*(gamm+H.*(4.*k.^2+l.^2).*eta)+tau)+H.*eta.*((k.^2+4.*l.^2).*U.*gamm+4.*H.*j2.^2.*U.*eta+j2.*tau));
t2=p.*ps;
trans=t1./t2;



return

end



