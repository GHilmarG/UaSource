

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
        iU=x<(lower+Eps) ; x(iU)=lower(iU)+Eps;
    else
        iU=x<(lower+Eps) ; x(iU)=lower+Eps;
    end
    
    if numel(upper)>1
        iL=x>(upper-Eps) ; x(iL)=upper(iL)-Eps;
    else
        iL=x>(upper-Eps) ; x(iL)=upper-Eps;
    end
    
end
