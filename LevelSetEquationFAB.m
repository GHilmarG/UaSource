function [k,dkdx]=LevelSetEquationFAB(CtrlVar,x,mu)
    
    

    % CtrlVar.LevelSetFABCostFunction="p2-q4" ;
    % mu=0.01*mu;
    % CtrlVar.LevelSetFABCostFunction="p4q4" ; 
    % CtrlVar.LevelSetFABCostFunction="p2q4"; 
    switch CtrlVar.LevelSetFABCostFunction
        
        case "p2q1"
            
            
            k=1-1./x ;
            k=k.*mu;
            
            if nargout>1
                dkdx=1./x.^2;
                dkdx=mu.*dkdx;
            end
            
        case "p2q2"
            
            xa=1;
            % P = 0.5*(x^2-xa^2)^2
            
            k=x.^2-xa^2 ;
            k=k.*mu;
            
            if nargout>1
                dkdx=2*x;
                dkdx=mu.*dkdx;
            end
            
        case "p2q4"
            
            xa=1; 
            xa=-CtrlVar.LSFslope;
            k=(x.^2-xa.^2).^3 ;
            k=k.*mu;
            
            if nargout>1
                dkdx=6*(x.^2-xa.^2).^2 .*x ;
                dkdx=mu.*dkdx;
            end
            
        case "p4q4"
            
            %if CtrlVar.time>100
            %   xa=1+CtrlVar.time/1000;
            %else
               xa=1;
               xa=-CtrlVar.LSFslope;
            %end
            
            k=(x.^4-xa.^4).^3 .* x.^2;
            k=k.*mu;
            
            if nargout>1
                dkdx=3.*(x.^4-xa.^4).^2 .*4.*x.^3 .* x.^2 + (x.^4-xa.^4).^3.*2.*x ;
                dkdx=mu.*dkdx;
            end
            
        case "Lietal2010"
            
            LT=x<1 ;
            k=1-1./x ;
            k(LT)=sin(2*pi*x(LT))./(2*pi*x(LT)+eps) ; % neg, then pos diffusion
            k=k.*mu;
            
            if nargout>1
                dkdx=1./x.^2;
                dkdx(LT)=2*pi*cos(2*pi*x(LT))./(2*pi*x(LT))-2*pi*sin(2*pi*x(LT))./(2*pi*x(LT)+eps).^2 ;
                dkdx=mu.*dkdx;
            end
            
            
            
    end
    
    
end