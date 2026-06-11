


function trans = T_SC_3vct(kx,ky,C,ca,t)

narginchk(5,5)

% Full Stokes, transient, basal slipperiness to surface topography transfer function 
% 

%% Original fortran code:
%
%    complex*16 function T_ZC_3uct(t,k,m,C,ca,Sh,Ch)
%       implicit none
%       complex*16 Im
% !      parameter (Im=cmplx(0.0d0,1.0d0)) g77 can't do it this way
% 
% c     History:
% c     first version using expressions from Maple written Nov-24-1994
% 
% 
%       real*8 m,k,C,ca,Sh,Ch,t
%       real*8 b,d,e
%       real*8 func_e,func_b,func_d,w_w,w_d
%       real*8 wd,ww
%       external func_a,func_b,func_d,w_d,w_w
% 
%       Im=cmplx(0.0d0,1.0d0)
% 
%       b=func_b(k,m,C,Sh,Ch,ca)
%       d=func_d(k,m,C,Sh,Ch)
%       e=func_e(k,m,C,Sh,Ch)
%       wd=w_d(k,m,C,ca,Sh,Ch)
%       ww=w_w(k,m,C,Sh,Ch)
% 
%       T_ZC_3uct=e/(d+Im*b)*(1.0d0-exp(Im*ww*t)*exp(-wd*t))
% c      T_ZC_3uct=e/(d+Im*b)*(1.0d0-exp(Im*d*t/cc)*exp(-b*t/cc))
% 
% c       T_ZC_3uct=C*T_ZC_3uct ! changed on May 5, 1997! and changed back on June 5 2000
%       return
%       end
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
%a=func_a(k,m,C,Sh,Ch);
b=func_b(k,m,C,Sh,Ch,ca);
d=func_d(k,m,C,Sh,Ch);
e = func_e(k,m,C,Sh,Ch) ; 
wd=w_d(k,m,ca,C,Sh,Ch);
ww=w_w(k,m,C,Sh,Ch);


T_ZC_3uct=(e./(d+Im*b)).*(1-exp(Im*ww.*t).*exp(-wd.*t)) ;

trans=T_ZC_3uct;


return
end

