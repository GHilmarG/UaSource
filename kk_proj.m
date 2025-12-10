

function [x,iU,iL] = kk_proj(x,upper,lower,Eps)
    %
    % projection onto active set
    %
    %ndim=length(x);
    %px=zeros(ndim,1);
    
    
    
    if nargin==3
        Eps=0;
    end
    
    if numel(lower) > 1
        iL=x<(lower+Eps) ; x(iL)=lower(iL)+Eps;
    else
        iL=x<(lower+Eps) ; x(iL)=lower+Eps;
    end
    
    if numel(upper)>1
        iU=x>(upper-Eps) ; x(iU)=upper(iU)-Eps;
    else
        iU=x>(upper-Eps) ; x(iU)=upper-Eps;
    end
    
end
