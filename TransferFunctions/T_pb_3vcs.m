function trans = T_pb_3vcs(kx,ky,C,ca)
    
    
    
    % Basal pressure, 3D full-Stokes, steady state
    
    if isvector(kx) && isvector(ky)
        % if k & l are vectors then
        k=repmat(kx,1,length(ky));
        l=repmat(ky',length(kx),1);
    else
        k=kx ; l=ky;
    end
    
    
    m=sqrt(k.^2+l.^2);
    t1 = 1i;
    t3 = C.^2;
    t6 = exp(4.*m);
    t7 = -1+t6;
    t8 = m.^2;
    t13 = exp(2.*m);
    t18 = t1.*t7;
    t23 = k.*C;
    t24 = 1+C;
    t25 = C.*t6;
    t27 = 2.*t13.*C;
    t28 = 4.*t13;
    t29 = -t25+t27+t28-C;
    t31 = t8.*m;
    t36 = t23.*t7.*t24.*t8;
    t48 = t24.^2;
    
    
    D = -t1.*k.*((t1.*t3.*t7.*t8+t1.*C.*(2.*t13+1+t6).*m+t18.*(2+C)).*ca-t23.*t24.*t29.*t31+t36-k.*t24.*(-C-t25+t27-2.*t6-2-t28).*m);
    N = (-t1.*t29.*m+t18).*ca+4.*k.*t13.*t48.*t31+t36+k.*(1+t6+6.*t13).*t24.*m;
    
    
    trans=D./N;
    
    
    
end
