

function fb = func_b(k,m,C,Sh,Ch,ca)


%%
% real*8 function func_b(k,m,C,Sh,Ch,ca)
%
%    implicit none
%    real*8 k,m,C,Sh,Ch,ca
%
%    func_b = ((m*C*Sh+Ch)*Sh-m)*ca
%
%
%    return
%    end
%

%%

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

fb = ((m.*C.*Sh+Ch).*Sh-m).*ca ;


%%


end