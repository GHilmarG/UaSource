function [detrended,a0,a1,a2,a3,a4]=detrend_xt(y,x,trend)
    
    a0=0 ; a1=0; a2=0; a3=0; a4=0;
    if nargin==2
        trend=1;
    end
    x=x(:) ; y=y(:);
    if trend==1
        sol=[ones(length(x),1) x]\y ;
        detrended=y-sol(1)-sol(2)*x;
        a0=sol(1); a1=sol(2);
    elseif trend==2
        sol=[ones(length(x),1) x x.*x]\y ;
        detrended=y-sol(1)-sol(2)*x-sol(3)*x.*x;
        a0=sol(1); a1=sol(2); a2=sol(3) ;
    elseif trend==3
        sol=[ones(length(x),1) x x.^2 x.^3]\y ;
        detrended=y-sol(1)-sol(2)*x-sol(3)*x.^2-sol(4)*x.^3;
        a0=sol(1); a1=sol(2); a2=sol(3) ; a3=sol(4);
    elseif trend==4
        sol=[ones(length(x),1) x x.^2 x.^3 x.^4]\y ;
        detrended=y-sol(1)-sol(2)*x-sol(3)*x.^2-sol(4)*x.^3-sol(5)*x.^4;
        a0=sol(1); a1=sol(2); a2=sol(3) ; a3=sol(4); a4=sol(5);
    end
    
    
    
    
end

