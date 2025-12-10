



function trans = T_SB_3vct(kx,ky,C,ca,t)

narginchk(5,5)

% Full Stokes, transient, bed topography to surface topography transfer function 

%% original fortran code 
% 
%   complex*16 function T_ZZ_3uct(t,k,m,C,ca,Sh,Ch)
% 
%       implicit none
%       complex*16 Im
%       parameter (Im=cmplx(0.0d0,1.0d0)) ! g77 can't do it this way
% 
%       real*8 m,k,C,ca,Sh,Ch,t
%       real*8 a,b,cc,d
%       real*8 func_a,func_b,func_d,w_w,w_d
%       real*8 wd,ww
% 
%       external func_a,func_b,func_d,w_w,w_d
% 
%       a=func_a(k,m,C,Sh,Ch)
%       b=func_b(k,m,C,Sh,Ch,ca)
%       d=func_d(k,m,C,Sh,Ch)
%       wd=w_d(k,m,C,ca,Sh,Ch)
%       ww=w_w(k,m,C,Sh,Ch)
%
%       T_ZZ_3uct=a/(d+Im*b)*(1.0d0-exp(Im*ww*t)*exp(-wd*t))
% c      T_ZZ_3uct=a/(d+Im*b)*(1.0d0-exp(Im*d*t/cc)*exp(-b*t/cc))
%
%       return
%       end
%
%%
if isvector(kx) && isvector(ky)
    % if k & l are vectors then
    k=repmat(kx,1,length(ky));
    l=repmat(ky',length(kx),1);
else
    k=kx ; l=ky;
end

m=sqrt(k.^2+l.^2);
Sh = sinh(m);
Ch = cosh(m);
Im=1i;
%%
a=func_a(k,m,C,Sh,Ch);
b=func_b(k,m,C,Sh,Ch,ca);
d=func_d(k,m,C,Sh,Ch);
wd=w_d(k,m,ca,C,Sh,Ch);
ww=w_w(k,m,C,Sh,Ch);

T_ZZ_3uct=(a./(d+Im*b)).*(1-exp(Im*ww*t).*exp(-wd*t));

%%

trans=T_ZZ_3uct;

return
end

