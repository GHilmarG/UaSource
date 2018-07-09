function [Tsb,Tsc,Tss0,Tub,Tuc,Tus0,Tvb,Tvc,Tvs0,Twb,Twc,Tws0] = ForwardFunctions2D(dx,dy,nx,ny,time,alpha,H,eta,C0,rho,g,m)
  
   
    k=fftspace(nx,dx); k(1)=eps ;
    l=fftspace(ny,dy); l(1)=eps ;
    Tss0=SSTREAM_Tss_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    Tsb=SSTREAM_Tsb_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    Tsc=SSTREAM_Tsc_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    
    Tus0=SSTREAM_Tus_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    Tub=SSTREAM_Tub_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    Tuc=SSTREAM_Tuc_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    
    
    Tvs0=SSTREAM_Tvs_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    Tvb=SSTREAM_Tvb_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    Tvc=SSTREAM_Tvc_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    
    
    Tws0=SSTREAM_Tss_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    Twb=SSTREAM_Tsb_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    Twc=SSTREAM_Tsc_t_3d_m(k,l,time,alpha,H,eta,C0,rho,g,m);
    
    
end

