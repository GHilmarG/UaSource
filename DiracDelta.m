function y = DiracDelta(k,x,x0)
	% Approximation to the Dirac Delta function
	%  The width of the approximation is about 1/k
	%  the limit k -> infty is the Dirac Delta  function
    
    if nargin==2 ; error('DiracDelta: Need three arguments \n') ; end
    
	
	y=2*k./(exp(k*(x-x0)) +exp(-k*(x-x0))).^2;
    
	
    y(isnan(y))=0; % underflow errors 
    
end

