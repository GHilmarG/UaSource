function [x,iActiveSet]=ProjGradient(x,dfdx,gamma,xmin,xmax)
    
    
	x=x+gamma*dfdx;
	x(x<xmin)=xmin;
    x(x>xmax)=xmax;
    
    
    iActiveSet= (x==xmin ) | (x==xmax)  ;
    
    
    %fprintf('Constraints active %-i, not active %-i \n ',sum(iActiveSet),numel(iActiveSet)-sum(iActiveSet))
end
    
    