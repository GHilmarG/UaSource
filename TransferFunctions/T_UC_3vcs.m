function trans = T_UC_3vcs(kx,ky,C,ca)

[k,l] = kxky2kl(kx,ky);

  t4=k.^2+l.^2;
  m=sqrt(t4);
 
      t1 = k.^2;
      t2 = C*t1;
%      t4 = m.^2;
      t7 = 1+C;
      t8 = t7*C;
      t9 = t4.^2;
      t11 = 2+C;
      t19 = sinh(m);
      t21 = cosh(m);
      t23 = t4*C;
      t34 = t19.^2;
      t37 = t7^2;
      t49 = C^2;
      t50 = t49*t4;
      t57 = t21.^2;
      t73 = t7*t4;
      t90 = C*(((t2-2).*t4*ca+1i*k.*(2*t8*t9-(t11*t7*t1+4).*t4+2*t1)).*t19.*t21+m.*((-2*t23+t1.*t11).*ca+1i*k.*((-t8*t1+2).*t4+t2)).*t34+2*m.*(t4*ca+1i*k.*(t37*t9+(-t37*t1+2+C).*t4-t1)))./t4./((1i*(t50+2).*k.*t7-3*C*ca).*m.*t57.*t21+((2+3*C)*ca-1i*((-2-2*C+t49)*t4-2).*k*t7).*m.*t21+((-2-t50)*ca+3*1i*C*k.*t73).*t34.*t19+((t23-2).*ca+1i*C*t7*(t73+4).*k.*t4).*t19);

  trans=t90;

return


