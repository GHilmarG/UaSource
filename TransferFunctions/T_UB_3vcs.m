function trans = T_UB_3vcs(kx,ky,C,ca)
    
    %  Transferfunction, multiplication with Zb(k,l) gives the
    %  x-component of the surface velocity field.
    %  3D, surface, constant viscosity, steady state
    
    
    %     History:
    %     first version using expressions from Maple written Nov-8-1994
    
    [k,l] = kxky2kl(kx,ky);
    
    
    t9=k.^2+l.^2;
    m=sqrt(t9);
    
    if max(m) > 220 ;
        disp(' T_UB_3vcs : m this large has been known in the past to cause numerical overflow problems ')
    end
    
    
    t1 = 1+C;
    t3 = k.^2;
    t7 = m.^2;
    t8 = t7.^2;
    t11 = C^2;
    t16 = t3.*t1;
    t20 = t11.*t7;
    t29 = sinh(m);
    t31 = cosh(m);
    t44 = t7.*m;
    t51 = t11.*t3;
    t52 = 2+C;
    t56 = t52.^2;
    t62 = t29.^2;
    t85 = t31.^2;
    t101 = t7.*t1;
    t118 = ((1i.*((3.*C.*t1.*t3-4).*k.*t8.*C+2.*(((4.*C+2+t11).*t3+4+4.*C).*t7-2.*t16).*k)+(((-t20-4-C).*t3+2).*C+4).*t7.*ca).*t29.*t31+(1i.*(((t16+2).*t11.*t7.*C+C.*(5.*C+4).*t3+2.*(3.*C-2).*t1).*k.*t44-4.*t3.*k.*m.*t1.*C)+((-3.*t51+2.*C.*t52).*t7-t3.*t56).*m*ca).*t62+2.*(-t51-C-2).*t44.*ca+2.*1i.*(-(2+5.*C+4.*t11).*t8+((5.*t11+6.*C+2).*t3-4-4.*C).*t7+2.*t16).*k.*m)./t7./((1i.*(t20+2).*k.*t1-3.*C.*ca).*m.*t85.*t31+((2+3.*C).*ca-1i.*((-2-2.*C+t11).*t7-2).*k.*t1).*m.*t31+((-2-t20).*ca+3.*1i.*C.*k.*t101).*t62.*t29+((t7.*C-2).*ca+1i.*(t101+4).*C*t1.*k.*t7).*t29);
    
    trans=t118;
    
    
end








