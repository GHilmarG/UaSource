function [xmin] = QuadtradicFit(Slope0,f0,l,fl)
    
    % finds the minimum based on a quadrdic fit given two function values at zero and l, and a slope at 0
    % 
    
	a2=(fl-f0-Slope0*l)/l^2;
	xmin=-Slope0/2/a2;
	
end

