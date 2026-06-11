function trans = T_VB_3vcs(kx,ky,C,ca)
  
  

if isvector(kx) & isvector(ky)
  % if k & l are vectors then 
  k=repmat(kx,1,length(ky));
  l=repmat(ky',length(kx),1);
else
  k=kx ; l=ky;
end

  t3=k.^2+l.^2;
  m=sqrt(t3);
  

  t2 = C^2;
%  t3 = m^2;
  t4 = t2*t3;
  t9 = 1+C;
  t11 = t3.^2;
  t18 = sinh(m);
  t20 = cosh(m);
  t23 = (2+C)^2;
  t38 = t18.^2;
  t59 = t20.^2;
  t75 = t3.*t9;
  
  t93 = k.*l.*(((-t4-4-C).*t3*C*ca+i*(3*t2*t9*t11+2*(4*C+2+t2)*t3-4-4*C).*k).*t18.*t20+((-3*t4-t23).*m*ca+i*(t2*C*t9*t11+C*(5*C+4)*t3-4*C*t9).*k.*m).*t38-2*t2*t3.*m.*ca+2*i*k.*m.*((C*(5*C+6)+2)*t3+2+2*C))./t3./((i*(t4+2).*k.*t9-3*C*ca).*m.*t59.*t20+((2+3*C)*ca-i*((-2-2*C+t2)*t3-2).*k.*t9).*m.*t20+((-2-t4)*ca+3*i*C*k.*t75).*t38.*t18+((t3*C-2).*ca+i*(t75+4).*C*t9.*k.*t3).*t18);

  
  trans=t93;
  
  return

