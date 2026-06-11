function trans = T_SC_3vcs(kx,ky,C,ca)
% kx and ky must be vectors


if isvector(kx) & isvector(ky)
  % if k & l are vectors then 
  k=repmat(kx,1,length(ky));
  l=repmat(ky',length(kx),1);
else
  k=kx ; l=ky;
end

t11=k.^2+l.^2;
m=sqrt(t11);

t2 = cosh(m);
t5 = 1.0+C;
t7 = sinh(m);
t9 = C*m.*t7+t2;
%t11 = m^2
t23 = -C*k.*m.*t2./(t5*m.*k.*(t2.*t9+1+t11.*t5)+i*(t9.*t7-m)*ca);
trans=t23;

return



