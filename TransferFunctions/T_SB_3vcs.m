function trans = T_SB_3vcs(kx,ky,C,ca)
    
    % FS, steady-state non-dimensional
    
    if isvector(kx) && isvector(ky)
        % if k & l are vectors then
        k=repmat(kx,1,length(ky));
        l=repmat(ky',length(kx),1);
    else
        k=kx ; l=ky;
    end
    
    t9=k.^2+l.^2;
    m=sqrt(t9);
    
    
    t1 = 1+C;
    t3 = sinh(m);
    t5 = cosh(m);
    t6 = C*m.*t3+t5;
    t8 = C^2;
    % t9 = m.^2
    t27 = ((t1.*t6+(1+C+t8.*t9).*t5).*m.*k)./(m.*k.*t1.*(t5.*t6+1+t9.*t1)+1i*(t6.*t3-m)*ca);
    trans=t27;
    
    
end

