function [R,K]=fRK(x,problemtype)

x=x(:) ;

% nargout

switch problemtype

    case "[x1,x2]"

        R=x  ;
        K=eye(numel(x)) ;


    case "[x1^2,x2^2]"

        R=x.^2;
        K=2*[x(1)  0     ;
            0   x(2) ] ;

    case "[x1^-100 x1]"

        R(1)=(x(1)^4-100*x(1)^2) ;
        R(2)=0;
        R=R(:) ;
        K=[4*x(1)-200*x(1) 0 ; ...
            0              0 ] ;

    case "[x1^-100 x1,x2^2]"

        R(1)=(x(1)^4-100*x(1)^2) ;
        R(2)=x(2)^2;
        R=R(:) ;
        K=[4*x(1)-200*x(1) 0 ; ...
            0              2*x(2) ] ;


    case "[x1+x2,x2]"

        R(1)=x(1)+x(2);
        R(2)=x(2);
        R=R(:) ;
        K=[1  1 ; 0 1] ;

    case "[x1^2,x2]"

        R(1)=x(1)^2;
        R(2)=x(2);
        R=R(:) ;
        K=[2*x(1)  0 ; 0 1] ;


    case "[x1^2+x2,x2]"

        R(1)=x(1)^2+x(2);
        R(2)=x(2);
        R=R(:) ;

        K=[2*x(1)  1 ; ...   % \nabla R1^T
            0     1] ;       % \nabla R2^T

    case "[x1^2+x2,x2^2+x1]"

        R(1)=x(1)^2+x(2);
        R(2)=x(2)^2+x(1);
        R=R(:) ;

        K=[2*x(1)   1   ; ...   % \nabla R1^T
            1   2*x(2)] ;       % \nabla R2^T

    case "[x1^3-100 x2,-x2^2+10 x1]"

        R(1)=x(1)^3-100*x(2);
        R(2)=-x(2)^2+10*x(1);
        R=R(:) ;

        K=[3*x(1)^2   -100   ; ...   % \nabla R1^T
            +10   -2*x(2)] ;         % \nabla R2^T

    case "Rosenbrock"

        [f,g,H] = RosenbrockFunction(x) ;

        R=g ; K = H ;

    case "lsqRosenbrock"

        [R,K] = lsqRosenbrockFunction(x) ;

    case "Beale"


        % R is 3
        R(1)= (1.5-x(1)+x(1).*x(2)).^2 ; 
        R(2)= (2.25-x(1)+x(1).*x(2).^2).^2 ; 
        R(3)=(2.625-x(1)+x(1).*x(2).^3).^2 ;
        R=R(:) ; 
        % K is 2 x 3
        K=[2*(1.5-x(1)+x(1).*x(2)).*(x(2)-1)            ,  2* (1.5-x(1)+x(1).*x(2))    .*x(1)                    ;   ...
           2*(2.25-x(1)+x(1).*x(2).^2).*(x(2).^2-1)     ,  2*(2.25-x(1)+x(1).*x(2).^2)  .*2*x(1).*x(2)           ;   ...
           2*(2.625-x(1)+x(1).*x(2).^3).*(x(2).^3-1)    ,  2*(2.625-x(1)+x(1).*x(2).^3).*3*x(1).*x(2).^2 ]  ;


    otherwise

        error("case not found")


end


R=R(:);

end

