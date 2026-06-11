function fd = func_d(k,m,C,Sh,Ch)

% real*8 function func_d(k,m,C,Sh,Ch)
% 
%       implicit none
%       real*8 k,m,C,Sh,Ch
%       real*8 t2,t9,t13
% 
%       t2 = 1.0d0+C
%       t9 = m**2
%       t13 = m*k*t2*(Ch*(m*C*Sh+Ch)+1.0d0+t9*t2)
% 
%       func_d=t13
%       return
%       end
% %%

% if isvector(kx) && isvector(ky)
%     % if k & l are vectors then
%     k=repmat(kx,1,length(ky));
%     l=repmat(ky',length(kx),1);
% else
%     k=kx ; l=ky;
% end
% 
% m=sqrt(k.^2+l.^2);
% Sh = sinh(m);
% Ch = cosh(m);


%%


t2 = 1+C;
t9 = m.^2;
t13 = m.*k.*t2.*(Ch.*(m.*C.*Sh+Ch)+1+t9.*t2);

fd=t13; 

%%


end