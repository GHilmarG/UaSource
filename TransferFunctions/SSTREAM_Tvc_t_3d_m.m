function [trans]=SSTREAM_Tvc_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m)


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
ps=m.*gamm.^2+H.*(l.^2.*(4+m)+k.^2.*(1+4.*m)).*gamm.*eta+4.*H.^2.*j2.^2.*eta.^2;
%xi=gamm+4.*H.*j2.*eta;
%nu=gamm+H.*(k.^2+4.*l.^2).*eta;
p=SSTREAM_Pole_t_3d_m(kx,ky,alpha,H,eta,C,rho,g,m);

if isnan(t)
  expp=0;
else
  expp=exp(p.*t);
end

%t1=H.*k.*l.*U.*gamm.*(-3i.*k.*U.*eta+ca.*tau-expp.*(3.*p.*eta-3i.*k.*U*eta+ca*tau));
%t2=p.*(gamm+H.*j2.*eta).*xi;
 
t1=H.*k.*l.*U.*gamm.*(-3i.*k.*U.*eta+ca.*tau-expp.*(3*(p.*eta-1i.*k.*U.*eta)+ca.*tau));
t2=p.*ps;


trans=t1./t2;


return

end



