



function fa = func_a(k,m,C,Sh,Ch)

narginchk(5,5)

%%  original fortran code
%
% real*8 function func_a(k,m,C,Sh,Ch)
%
%    implicit none
%    real*8 k,m,C,Sh,Ch
%
%    real*8 t8,t9,t15
%
%    t8 = C**2
%    t9 = m**2
%    t15 = ((1.0d0+C)*(m*C*Sh+Ch)+(1.0d0+C+t8*t9)*Ch)*m*k
%
%    func_a=t15
%
%    return
%    end
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

t8 = C.^2;
t9 = m.^2;
t15 = ((1+C).*(m.*C.*Sh+Ch)+(1+C+t8.*t9).*Ch).*m.*k;

fa=t15;


%%


end