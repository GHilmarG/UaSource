function [R,K] = lsqRosenbrockFunction(x)


% minimum is at (1,1)

% Calculate objective f

c=sqrt(100);
R(1) = c* (x(2) - x(1).^2)  ;
R(2) =  1-x(1) ;

% f = (sqrt(100)*(x(2) - x(1).^2)).^2 + (1-x(1)).^2;

if nargout > 1 % gradient required

    K=[-2*c*x(1)  ,    c ; ...
        -1        ,     0 ] ;


end



end