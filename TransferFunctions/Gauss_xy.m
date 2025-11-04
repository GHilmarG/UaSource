




function [s,u,v,w,db,dc]=Gauss_xy(x,dx,y,dy,alpha,C,ampl_b,sigma_bx,sigma_by,ampl_c,sigma_cx,sigma_cy,theta,t)

% Calculates surface shape (s) and velocity (u,w) as linear medium
% flows over a Gaussian shaped bedrock and slipperiness perturbations.

% 3d map-pane (surface fields only)

% alpha is mean slope, C is mean slipperiness, 
% x/y and dx/dy model range and grid spacing,
% x and y are vectors 
% ampl_b/c and sigma_b/c are the amplitudes and the widths of the Gaussian perturbations

% the shape of the Gaussian peak is determined by the parameters
% amp_b,sigma_bx, sigma_by and theta where:
% amp is the amplitude, sigma_x and sigma_y the widths, and theta the
% orientation of the peak with respect to mean flow direction

% _b refers to a B pert, and _c to a C pert


if nargin == 13

    t=nan ;
    transient=false;
else
    transient=true;
end


[X,Y]=meshgrid(x,y);
nx=length(x) ; ny=length(y) ; ca=cot(alpha);

Xt=X*cos(theta)+Y*sin(theta);
Yt=-X*sin(theta)+Y*cos(theta);

% create bed and slipperiness perturbations and do forward fft
db=ampl_b*exp(-Xt.^2./sigma_bx^2-Yt.^2./sigma_by^2); db=db-mean(db(:)) ; Delta_b=ifft2(db);
dc=ampl_c*exp(-Xt.^2./sigma_cx^2-Yt.^2./sigma_cy^2); dc=dc-mean(dc(:)) ; Delta_c=ifft2(dc);


kx=fftspace(nx,dx); kx(1)=eps ;
ky=fftspace(ny,dy); ky(1)=eps ;

% multiply transfer functions with corresponding basal perturbations
% do inverse fft and only keep the real part

if transient


    s=real(fft2(T_SB_3vct(kx,ky,C,ca,t).*Delta_b+T_SC_3vct(kx,ky,C,ca,t).*Delta_c));
   
    % for the transient case, velocity transfer functions have not been implemented
    %u=real(fft2(T_UB_3vct(kx,ky,C,ca,t).*Delta_b+T_UC_3vct(kx,ky,C,ca,t).*Delta_c));
    %v=real(fft2(T_VB_3vct(kx,ky,C,ca,t).*Delta_b+T_VC_3vct(kx,ky,C,ca,t).*Delta_c));
    %w=real(fft2(T_WB_3vct(kx,ky,C,ca,t).*Delta_b+T_WC_3vct(kx,ky,C,ca,t).*Delta_c));
    u=nan ; v=nan ; w=nan ;

else

    s=real(fft2(T_SB_3vcs(kx,ky,C,ca).*Delta_b+T_SC_3vcs(kx,ky,C,ca).*Delta_c));
    u=real(fft2(T_UB_3vcs(kx,ky,C,ca).*Delta_b+T_UC_3vcs(kx,ky,C,ca).*Delta_c));
    v=real(fft2(T_VB_3vcs(kx,ky,C,ca).*Delta_b+T_VC_3vcs(kx,ky,C,ca).*Delta_c));
    w=real(fft2(T_WB_3vcs(kx,ky,C,ca).*Delta_b+T_WC_3vcs(kx,ky,C,ca).*Delta_c));

end


return

end


