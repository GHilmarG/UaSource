function B=UaResize(A,nB,varargin)

%%
% 
%  returns a vector with length nB
%
%  if on input the vector A has fewer elements than nB, A is padded with nan 
%
%%


nA=length(A);


if nA >= nB

    B=A ;
    

else

    B=nan(nB,1);
    B(1:nA)=A(:) ; 



end






end
