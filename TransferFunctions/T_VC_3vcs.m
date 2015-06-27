function trans = T_VC_3vcs(kx,ky,C,ca)


if isvector(kx) & isvector(ky)
  % if k & l are vectors then 
  k=repmat(kx,1,length(ky));
  l=repmat(ky',length(kx),1);
else
  k=kx ; l=ky;
end

t3=k.^2+l.^2;
m=sqrt(t3);


t1 = C*k;
%      t3 = m^2
t4 = t3*C;
t7 = 1+C;
t13 = sinh(m);
t15 = cosh(m);
t20 = t3*t7;
t26 = t13.^2;
t28 = t7^2;
t36 = C^2;
t37 = t36*t3;
t44 = t15.^2;
t76 = -t1.*l.*((-t4*ca+i*((2+C)*t7*t3-2).*k).*t13.*t15+((-2-C)*m*ca+i*(t20-1).*k.*m*C).*t26+2*i*(t3.*t28+1).*m.*k)./t3./((i*(t37+2).*k.*t7-3*C*ca).*m.*t44.*t15+((2+3*C)*ca-i*((-2-2*C+t36)*t3-2).*k*t7).*m.*t15+((-2-t37)*ca+3*i*t1.*t20).*t26.*t13+((t4-2)*ca+i*C*t7*(t20+4).*k.*t3).*t13);

trans=t76;
