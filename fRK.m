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


    otherwise

        error("case not found")


end