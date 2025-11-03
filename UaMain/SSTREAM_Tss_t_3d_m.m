function [trans]=SSTREAM_Tss_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m)
    
    
    % time-dependent ratio between surface topography and
    % bedrock in Fourier space in dimensional units
    
    % for non-dimensional a la gudmundsson 2003
    % put eta=1/2, H0=1,rho*g=1/sin(alpha), lambda in units of mean ice thickness
    % Cnondimensional=2 eta/H0 Cdimesional, and then
    % U0nondimensional=Cnondimensional
    
    
    
    p=SSTREAM_Pole_t_3d_m(kx,ky,alpha,H,eta,C,rho,g,m);
    
    if isnan(t)
        expp=0;
        %disp(' steady state option used ')
    else
        expp=exp(p*t);
    end
    
    trans=expp;
    
    
    return
    
end



