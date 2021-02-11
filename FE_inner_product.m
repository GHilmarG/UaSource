function [fg,theta]=FE_inner_product(f,g,M)



nf=numel(f) ; ng=numel(g) ;

[nM,mM]=size(M);

if nM==nf
    
    fg=f'*M*g;   % = (M*f)' * g
    fNorm=sqrt(f'*M*f) ; 
    gNorm=sqrt(g'*M*g) ; 
    theta=acos(fg/(fNorm*gNorm)) ; 
    
elseif nf==2*nM
    
    fg=f(1:nM)'*M*g(1:nM)+f(1+nM:nf)'*M*g(1+nM:ng);
    
    fNorm=sqrt(f(1:nM)'*M*f(1:nM)+f(1+nM:nf)'*M*f(1+nM:nf)) ;
    gNorm=sqrt(g(1:nM)'*M*g(1:nM)+g(1+nM:ng)'*M*g(1+nM:ng)) ;
    theta=acos(fg/(fNorm*gNorm)) ; 
    
end




end



