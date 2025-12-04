
function fc = func_c(k,m,C,Sh,Ch)

%%
%    real*8 function func_c(k,m,C,Sh,Ch)
% 
%       implicit none
%       real*8 k,m,C,Sh,Ch
% 
%       real*8 t1,t7,t10,t12
% 
% c     wrong, error discovered on August 15, 1995
% c      t1 = m**2
% c      t10 = (t1*(1.0d0+C)+m*C*Sh+Ch)*m*Ch
% 
% 
%       t1 = m**2
% c      t7 = cosh(m)
%       t7 = Ch
%       t10 = t7**2
% c      t12 = t1*m*(C+1)+t1*C*sinh(m)*t7+m*t10
%       t12 = t1*m*(C+1)+t1*C*Sh*t7+m*t10
% 
%       func_c=t12
%       return
%       end
%%

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

t1 = m.^2; 
t10 = Ch.^2 ;
t12 = t1.*m.*(C+1)+t1.*C.*Sh.*Ch+m.*t10 ;

fc = t12 ; 




%%


end