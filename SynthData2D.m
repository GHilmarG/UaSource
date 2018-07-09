function [s,u,v,w,db,dc,ds0] = SynthData2D(x,y,time)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
hmean=1000;   
ampl_b=0.1*hmean     ; sigma_bx=5*hmean ; sigma_by=5*hmean ;
ampl_c=0.1           ; sigma_cx=5*hmean ; sigma_cy=5*hmean ;
ampl_s=0.00          ; sigma_sx=5*hmean ; sigma_sy=5*hmean ;

rho=917; g=9.81/1000; alpha=0.01;
dt=0;

m=1 ; n=1 ;
SlipRatio=1000; ud=1 ;

taub=rho*g*hmean*sin(alpha);

%AGlen=1.0e-15*1e3*365.25*24*60*60;
%ud=2*AGlen*taub^n*hmean/(n+1);

AGlen=ud/(2*taub^n*hmean/(n+1));
ub=SlipRatio*ud;

C0=ub/taub^m;

[X,Y]=ndgrid(x,y);

nx=length(x) ; ny=length(y) ; dx=x(2)-x(1);dy=y(2)-y(1);
theta=0 ;

Xt=X*cos(theta)+Y*sin(theta);
Yt=-X*sin(theta)+Y*cos(theta);


ds0=ampl_s*exp(-Xt.^2./sigma_sx^2-Yt.^2./sigma_sy^2); ds0=ds0-mean(ds0(:)) ; delta_s0=ifft2(ds0);
db=ampl_b*exp(-Xt.^2./sigma_bx^2-Yt.^2./sigma_by^2); db=db-mean(db(:)) ; delta_b=ifft2(db);
dc=ampl_c*exp(-Xt.^2./sigma_cx^2-Yt.^2./sigma_cy^2); dc=dc-mean(dc(:)) ; delta_c=ifft2(dc);


k=fftspace(nx,dx); k(1)=eps ;
l=fftspace(ny,dy); l(1)=eps ;

H=hmean; eta=1/(2*AGlen);

s=real(fft2(SSTREAM_Tss_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m).*delta_s0)+...
    fft2(SSTREAM_Tsb_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m).*delta_b)+...
    fft2(SSTREAM_Tsc_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m).*delta_c));

u=real(fft2(SSTREAM_Tus_t_3d_m(k,l,time-dt,alpha,H,eta,C0,rho,g,m).*delta_s0)+...
    fft2(SSTREAM_Tub_t_3d_m(k,l,time-dt,alpha,H,eta,C0,rho,g,m).*delta_b)+...
    fft2(SSTREAM_Tuc_t_3d_m(k,l,time-dt,alpha,H,eta,C0,rho,g,m).*delta_c));


v=real(fft2(SSTREAM_Tvs_t_3d_m(k,l,time-dt,alpha,H,eta,C0,rho,g,m).*delta_s0)+...
    fft2(SSTREAM_Tvb_t_3d_m(k,l,time-dt,alpha,H,eta,C0,rho,g,m).*delta_b)+...
    fft2(SSTREAM_Tvc_t_3d_m(k,l,time-dt,alpha,H,eta,C0,rho,g,m).*delta_c));


w=real(fft2(SSTREAM_Tws_t_3d_m(k,l,time-dt,alpha,H,eta,C0,rho,g,m).*delta_s0)+...
    fft2(SSTREAM_Twb_t_3d_m(k,l,time-dt,alpha,H,eta,C0,rho,g,m).*delta_b)+...
    fft2(SSTREAM_Twc_t_3d_m(k,l,time-dt,alpha,H,eta,C0,rho,g,m).*delta_c));

s=s+hmean ;
u=u+ub ;




    
end

