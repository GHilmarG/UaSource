function y = TopHatApprox(k,x,halfwidth)
    
    % Creates a smooth top-hat function centered around x=0 with width of about 2*halfwidth
    % and a slopes of approx widths of 1/k 
    %
    % x=linspace(-100,100,1000) ; y = TopHatApprox(10,x,20) ; figure ; plot(x,y)
    % x=linspace(-100,100,1000) ; y = TopHatApprox(1/10,x,20) ; figure ; plot(x,y)
    
    y = 1 - HeavisideApprox(k,x,halfwidth) - HeavisideApprox(k,-x,halfwidth) ;
    
    
    
    
    
    
end