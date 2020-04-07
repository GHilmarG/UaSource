function  [r,rbcs,ruv,rh]=ResidualCostFunction(CtrlVar,MUA,L,Rfields,Rbcs,fext0,uvORuvh)
                          


nargoutchk(1,4)
narginchk(7,7)

% square of normalized residual
% returns (R1'*R1+R2'*R2)/(F0'*F0);
%
% ruv and rh are normalized with F0'F0
% rl is not normalized

if isempty(Rbcs) ; Rbcs=0 ; end

rfields=Rfields'*Rfields;
rbcs=Rbcs'*Rbcs;

% TestIng
fext0=fext0+eps;

if contains(uvORuvh,"-uvh-")
    
   % the h part of Fext0 goes to zero with dt, so redefine and make independent of dt 
    fext0(2*MUA.Nnodes+1:end)=fext0(2*MUA.Nnodes+1:end)/CtrlVar.dt ;
    
end



% r=full(real((rfields+rbcs)/(Fext0'*Fext0)));
r=full(real(rfields/(fext0'*fext0)));  % August 2017. Decided to include only the uvh fields in the residual. 
                                       % After all the BCs are linear and are always fullfilled exactly.
                                       % The rbcs residulas will always be
                                       % tiny (e.g 1e-50) excep at the
                                       % beginning at the nl-iteration
                                       % where the system has not been
                                       % solved. 


if nargout > 2
    
    Nuv=2*MUA.Nnodes;
    
    ruv=full(Rfields(1:Nuv)'*Rfields(1:Nuv)/(fext0(1:Nuv)'*fext0(1:Nuv)));
    ruv=real(ruv);
    
    if nargout>3
        rh=full(Rfields(Nuv+1:end)'*Rfields(Nuv+1:end)/(fext0(Nuv+1:end)'*fext0(Nuv+1:end)));
        rh=real(rh);
    end
    
    
end



end

