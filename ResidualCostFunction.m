function  [r,rbcs,ruv,rh]=ResidualCostFunction(Rfields,Rbcs,Fext0,Nnodes)

nargoutchk(1,4)
narginchk(4,4)

% square of normalized residual
% returns (R1'*R1+R2'*R2)/(F0'*F0);
%
% ruv and rh are normalized with F0'F0
% rl is not normalized

if isempty(Rbcs) ; Rbcs=0 ; end

rfields=Rfields'*Rfields;
rbcs=Rbcs'*Rbcs;

%r=full(real((rfields+rbcs)/(Fext0'*Fext0)));
r=full(real(rfields/(Fext0'*Fext0)));  % August 2017. Decided to include only the uvh fields in the residual. 
                                       % After all the BCs are linear and are always fullfilled exactly.
                                       % The rbcs residulas will always be
                                       % tiny (e.g 1e-50) excep at the
                                       % beginning at the nl-iteration
                                       % where the system has not been
                                       % solved. 


if nargout > 2
    
    Nuv=2*Nnodes;
    
    ruv=full(Rfields(1:Nuv)'*Rfields(1:Nuv)/(Fext0(1:Nuv)'*Fext0(1:Nuv)));
    ruv=real(ruv);
    
    if nargout>3
        rh=full(Rfields(Nuv+1:end)'*Rfields(Nuv+1:end)/(Fext0(Nuv+1:end)'*Fext0(Nuv+1:end)));
        rh=real(rh);
    end
    
    
end



end

