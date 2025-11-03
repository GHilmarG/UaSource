function [f,g,H] = RosenbrockFunction(x)


% minimum is at (1,1)

% Calculate objective f

f1 = (sqrt(100)*(x(2) - x(1).^2)).^2 ;
f2 =  (1-x(1)).^2;

f=f1.^2+f2.^2 ; 

% f = (sqrt(100)*(x(2) - x(1).^2)).^2 + (1-x(1)).^2;

if nargout > 1 % gradient required
  
    g = [-400*(x(2)-x(1).^2).*x(1) - 2*(1-x(1));
        200*(x(2)-x(1).^2)];

    if nargout> 2
    
        H=zeros(2,2);

        H(1,1)=1200 * x(1).^2-400*x(2)+ 2; 
        H(1,2)=-400*x(1);
        H(2,1)=-400*x(1);
        H(2,2)= 200 ;

    end

end



end