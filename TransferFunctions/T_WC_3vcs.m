function trans = T_WC_3vcs(kx,ky,C,ca)


if isvector(kx) & isvector(ky)
  % if k & l are vectors then 
  k=repmat(kx,1,length(ky));
  l=repmat(ky',length(kx),1);
else
  k=kx ; l=ky;
end

t13=k.^2+l.^2;
m=sqrt(t13);

%     3D, surface, constant viscosity, steady state
%     This transferfunction uses the U_b = C (1 + dC) definition.

t1 = k.^2;
t2 = 1+C;
t5 = cosh(m);
t9 = sinh(m);
t11 = C*m.*t9+t5;
%      t13 = m.^2
t25 = i*t2*C*m.*t1.*t5./(t2*m.*k.*(t11.*t5+1+t13*t2)+i*(t11.*t9-m)*ca);

trans=t25;

return

