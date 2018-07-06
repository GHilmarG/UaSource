function y = HeavisideApprox(k,x,x0)
	
    if nargin==2 ; error('HeavisideApprox: Need three arguments \n') ; end
    
	% Approximation to the Heaviside step  function
	%  The width of the step is about 1/k
	%  the limit k -> infty is Heaviside step function
	%  
	y=1./(1+exp(-2*k*(x-x0)));
	
    
	
end