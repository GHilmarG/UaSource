function trans = T_WB_3vcs(kx,ky,C,ca)
  
    

if isvector(kx) & isvector(ky)
  % if k & l are vectors then 
  k=repmat(kx,1,length(ky));
  l=repmat(ky',length(kx),1);
else
  k=kx ; l=ky;
end

t11=k.^2+l.^2;
m=sqrt(t11);
t1 = k.^2;
t2 = 1+C;
t5 = sinh(m);
t7 = cosh(m);
t8 = C*m.*t5+t7;
t10 = C^2;
%  t11 = m.^2;
t30 = -i*t2*t1.*(t2*t8+(1+C+t10*t11).*t7).*m./(t2*m.*k.*(t7.*t8+1+t11*t2)+i*(t8.*t5-m)*ca);
  
  
trans=t30;
  
return

