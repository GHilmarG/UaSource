function [trans]=SSTREAM_Twc_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m)


% for non-dimensional a la gudmundsson 2003
% put eta=1/2, H0=1,rho.*g=1/sin(alpha), lambda in units of mean ice thickness
% Cnondimensional=2 eta/H0 Cdimesional, and then
% U0nondimensional=Cnondimensional 

[k,l] = kxky2kl(kx,ky);


j2=k.^2+l.^2;

tau=rho.*g.*sin(alpha).*H;
ca=cot(alpha);
U=C.*tau^m;
gamm=tau^(1-m) / (C.*m); % can't use gamma because it is a function
xi=gamm+4.*H.*j2.*eta;
nu=gamm+H.*(k.^2+4.*l.^2).*eta;
p=SSTREAM_Pole_t_3d_m(kx,ky,alpha,H,eta,C,rho,g,m);

if isnan(t)
  expp=0;
  %disp(' steady state option used ')
else
  expp=exp(p.*t);
end


% I have not yet tested this against numerical code
t1=H.*k.*U.*gamm.*(expp.*(1i.*p+k.*U)-k.*U);
t2=p.*xi;
trans=t1./t2;



return

end



