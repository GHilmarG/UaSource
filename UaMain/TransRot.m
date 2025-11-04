% trying to formulat the tranfer stuff in arbitrary coordinte system

% U W : velocities in spatial space and unrotated coordinate system
%
%   u = T_ub b
%
% [ U' ; W'] =  R(alpha) [U ; W]
%
%  [ u; w] = [F 0 ; 0 F] [U W]
%

%%

% Note : not finished
% I'm thinking it might be better just to solve this approximatly
% let dh be the thickness perturbation and define dh=ds+db with ds deviation from mean surface mesured
% in true veritcal direction, don't make a distinction beweeen horizontal velocity and veritcal in tilted and
% untitled coordinate system, just worry about the rotation in the xy plane


% Units: m kPa a
alpha=0.1*pi/180;

H=1000/tan(alpha) ;  % ice thickness measured normal to the mean surface
rho=917 ; g=9.81/1000;
Us=300 ; % m/a
n=1;
Ud=1;
m=3;
t=NaN;

tau=rho*g*H*sin(alpha); % kPa
%A=1.67e-7 ; % kPa^{-3} a^{-1} 0 degrees C
%Ud=2*A*tau^n *H /(n+1) ; % rough estimate for deformational velocity (m/a)

A=Ud*(n+1)/(2*H*tau^n);
Ub=Us-Ud;
SlipRatio=Ub/Ud;

eta=tau^(1-n)/2/A;  % effective viscosity

%Ub=C tau^m
C=Ub/tau^m;


nx=128; xmin=0 ; xmax=100;

dx=H;

x=[-dx*nx/2:dx:dx*(nx/2-1)]';  % model domains



dbAmpl=H/100; xb=-10 ; bsigma=10*H;
dcAmpl=0.1*C; xc=10 ; csigma=10*H;


b=-x*tan(alpha)+dbAmpl*exp( -(x-xb).^2/bsigma^2);
c=dcAmpl*exp( -(x-xc).^2/csigma^2);



R=sparse(1:2*nx,1:2*nx,cos(alpha),2*nx,2*nx);
R=R+sparse(1:nx,nx+1:2*nx,-sin(alpha),2*nx,2*nx);
R=R+sparse(nx+1:2*nx,1:nx,sin(alpha),2*nx,2*nx);


[temp] = R* [x;b];
xRotb=temp(1:nx) ; bRot=temp(nx+1:2*nx);



bbar=bRot*0+mean(bRot) ; db=bRot-bbar;



dsRot=real(fft(Fsb*ifft(db)));
duRot=real(fft(Fub*ifft(db)));
dwRot=real(fft(Fwb*ifft(db)));


kx=fftspace(nx,dx); kx(1)=eps ; % freq vector
lambda=2*pi./kx;

dsRot=real(fft(SSTREAM_Tsb_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*db+SSTREAM_Tsc_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*dc));
  duRot=real(fft(SSTREAM_Tub_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*db+SSTREAM_Tuc_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*dc));
  dwRot=real(fft(SSTREAM_Twb_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*db+SSTREAM_Twc_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*db));
%

% R' is the inverse of R

[test] = R'*[duRot ; dwRot];
du=test(1:nx); dw=test(nx+1:2*nx);

test=R'*[xRotb;dsRot];
ds=test(nx+1:2*nx); 
xback=test(1:nx) ; % this is not exacly equal to the original x, but the difference will be very small for small alpha

% in general I need to to a reinterpolation from xback to x, but even if I have the 2D
% problem this should not be expensive, and most certainly not really needed

figure ; plot(du) ; hold on ; plot(dw,'r')

%%
%  
%  db=ampl_b*exp(-(x-x_b).^2/sigma_b^2); db=db-mean(db) ; Delta_b=ifft(db);
%  dc=ampl_c*exp(-(x-x_c).^2/sigma_c^2); dc=dc-mean(dc) ; Delta_c=ifft(dc);
%  s=real(fft(SSTREAM_Tsb_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*Delta_b+SSTREAM_Tsc_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*Delta_c));
%  u=real(fft(SSTREAM_Tub_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*Delta_b+SSTREAM_Tuc_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*Delta_c));
%  w=real(fft(SSTREAM_Twb_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*Delta_b+SSTREAM_Twc_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m).*Delta_c));
%
 
 













        