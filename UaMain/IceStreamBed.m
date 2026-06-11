
function B=IceStreamBed(x,y,w,d,sigma)

    
    % w : half-width of ice stream bed
    % d : depth of the channel
    % sigma : sharpness of the channel (affects transvese slope)
    
    

B=720-778.5*x/750e3+d*(1./(1+exp(-2*(y-w)/sigma))+1./(1+exp(2*(y+w)/sigma)));


% B=d*(1./(1+exp(-2*(y-w)/sigma))+1./(1+exp(2*(y+w)/sigma)));
% dBdy=-2*d*exp(-2*(y-w)/sigma)./(1+exp(-2*(y-w)/sigma)).^2/sigma+...
%       2*d*exp(+2*(y+w)/sigma)./(1+exp(+2*(y+w)/sigma)).^2/sigma;
  
end



  