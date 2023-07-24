function [R,K]=fRK(x,problemtype)

x=x(:) ;

switch problemtype

    case "[x1,x2]"

        R=x + 1 ;
        K=eye(numel(x)) ;

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


    otherwise

        error("case not found")


end