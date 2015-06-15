function [pole]=SSTREAM_Pole_t_3d_m(kx,ky,alpha,H,eta,C,rho,g,m)
    
    [k,l] = kxky2kl(kx,ky);
        
    j2=k.^2+l.^2;
    ca=cot(alpha);
    tau=rho*g*sin(alpha)*H;
    U=C*tau^m;
    gamm=tau^(1-m) / (C*m); % can't use gamma because it is a function
    kappa=l.^2.*(4+m)+k.^2.*(1+4*m);
%    xi=gamm+4*H*j2*eta;   old, incorrect verstion used in TC paper
%    t1=(k.*(ca.*H.*k-1i)+ca.*H.*l.^2).*tau;
%    t2=xi;

% corrected expression	
    t1=tau*(-1i*k.*(m*gamm+H*j2*eta)+ca*H*((l.^2+k.^2*m)*gamm+H*j2.^2*eta));
    t2=m*gamm^2+H*eta*(4*H*j2.^2*eta+gamm*kappa);

    pole=1i*k*U-t1./t2;
    
    
    
    return
    
    
end
