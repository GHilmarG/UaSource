function ww = w_w(k,m,C,Sh,Ch)

narginchk(5,5)

%%
% 
%      real*8 function w_w(k,m,C,Sh,Ch)
% 
%       implicit none
%       real*8 k,m,C,Sh,Ch
% 
%       real*8 t1,t3,t5,t9,t12,t17
% 
%       t1 = C+1
%       t3 = m**2
%       t5 = t3*m*t1
% c      t9 = dcosh(m)
%       t9 = Ch
%       t12 = (m*C*Sh+t9)*m*t9
% c      t12 = (m*C*dsinh(m)+t9)*m*t9
%       t17 = k*t1*(t5+t12+m)/(t5+t12)
% 
%       w_w=t17
% 
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

t1 = C+1;
t3 = m.^2;
t5 = t3.*m.*t1;
t9 = Ch ; 
t12 = (m.*C.*Sh+t9).*m.*t9;
t17 = k.*t1.*(t5+t12+m)./(t5+t12);

ww=t17;

return
end

