function wd = w_d(k,m,ca,C,Sh,Ch)

narginchk(6,6)

%%
%
%    real*8 function w_d(k,m,C,ca,Sh,Ch)
%
%       implicit none
%       real*8 k,m,C,Sh,Ch,ca
%
%       real*8 t1,t2,t4,t7,t11,t20
%
%
%       t1 = Ch
% c      t1 = dcosh(m)
%       t2 = m*C
% c      t4 = dsinh(m)
%       t4 = Sh
%       t7 = C+1
%       t11 = m**2
%       t20 = (t1*(t2*t1+t4)-m*t7)*ca/(t11*m*t7+(t2*t4+t1)*m*t1)
%
%       w_d=t20
%       return
%       end
%%
% 
% 
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


t1 = Ch;
t2 = m.*C;
t4 = Sh;
t7 = C+1;
t11 = m.^2;
t20 = (t1.*(t2.*t1+t4)-m.*t7)*ca./(t11.*m.*t7+(t2.*t4+t1).*m.*t1) ;

wd=t20;


return
end

