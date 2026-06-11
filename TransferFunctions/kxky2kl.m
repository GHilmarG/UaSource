
function [k,l] = kxky2kl(kx,ky)
    % kxky2kl
    
    % if both kx and ky are vectors create 2d arrays by replicating
    % otherwise do nothing
    
    if isvector(kx) && isvector(ky) && length(kx) >1 && length(ky)>1
        kx=kx(:); ky=ky(:);
        k=repmat(kx,1,length(ky));
        l=repmat(ky',length(kx),1);
    else
        k=kx ; l=ky;
    end
    
end

