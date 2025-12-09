function y = DiracDelta(k,x,x0)

%%
% Approximation to the Dirac Delta function
%  The width of the approximation is about 1/k
%  the limit k -> infty is the Dirac Delta  function
%
%
%  To create a delta peak of width W and amplitude A centered around x0=0;
%
%  W=100 ; A=5; x0=0 ; x=linspace(-10*W,10*W,1000) ; Peak = A*2*W*DiracDelta(1/W,x,x0) ; figure ; plot(x,Peak)
%
% 
%
%  
% 
%  
%%

if nargin==2 ; error('DiracDelta: Need three arguments \n') ; end


y=2*k./(exp(k*(x-x0)) +exp(-k*(x-x0))).^2;


y(isnan(y))=0; % underflow errors

end

