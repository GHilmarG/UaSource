function dnewt=bfgsw(sstore,alpha,beta,ns,dsd,hess0)
    %
    % bfgsw
    %
    % C. T. Kelley, Dec 20, 1996
    %
    % This code comes with no guarantee or warranty of any kind.
    %
    % This code is used in bfgswopt.m
    %
    % There is no reason to ever call this directly.
    %
    % form the product of the bfgs approximate inverse Hessian
    % with a vector using the Steve Wright method
    %
    
    % dsd : negative of gradient
    
    userhess=0; 
    if nargin==6 ;
        userhess=1; 
    end
    
    dnewt=dsd;
    
    if userhess==1 ;
        dnewt=feval(hess0,dnewt); 
    end
    
    if (ns<=1) ; 
        return; 
    end
    
    dnewt=dsd; n=length(dsd);
    
    if userhess==1 ; 
        dnewt=feval(hess0,dnewt); 
    end
    
    sigma=sstore(:,1:ns-1)'*dsd; 
    gamma1=alpha(1:ns-1).*sigma;
    gamma2=beta(1:ns-1).*sigma;
    gamma3=gamma1+beta(1:ns-1).*(sstore(:,2:ns)'*dsd);
    delta=gamma2(1:ns-2)+gamma3(2:ns-1);
    dnewt=dnewt+gamma3(1)*sstore(:,1)+gamma2(ns-1)*sstore(:,ns);
    
    if(ns <=2) ; 
        return; 
    end
    
    dnewt=dnewt+sstore(1:n,2:ns-1)*delta(1:ns-2);
    
end
    
    %