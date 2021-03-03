function [fg,theta]=FE_inner_product(f,g,M)

%%
%
%   [fg,theta]=FE_inner_product(f,g,M)
%
%       f :     Nodal variables
%       g :     Nodal variables
%       M :     The FE mass matrix, this you can find as the field MUA.M 
%
%
%      fg :     FE inner product between the nodal variables f and g
%   theta :     The angle between f and g 
%
% Examples:
%
%   fg=FE_inner_product(f,g,MUA.M)   ; 
%
%   fNorm=sqrt(FE_inner_product(f,f,MUA.M))   ; 
% 
%
% Note:  If you just need to calculate the norm of a nodal variable f, you can do so simply as:
%
%   fNorm= f'*MUA.M* f 
%
% So you don't really need to call this m-File. 
%
%%


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



