function  [r,ruv,rh]=ResidualCostFunction(R,F0)
    
    % square of normalized residual
    % returns R'*R/(F0'*F0);
    
    if nargin==1
        r=full(R'*R); r=real(r) ; 
         if nargout > 1
            Nuv=2*numel(R)/3;
            ruv=full(R(1:Nuv)'*R(1:Nuv)); ruv=real(ruv); 
            rh=full(R(Nuv:end)'*R(Nuv:end)); rh=real(rh); 
        end
        
    else
        
        r=full(R'*R/(F0'*F0)); r=real(r) ; 
        
        if nargout > 1
            Nuv=2*numel(R)/3;
            ruv=full(R(1:Nuv)'*R(1:Nuv)/(F0(1:Nuv)'*F0(1:Nuv))); ruv=real(ruv); 
            rh=full(R(Nuv:end)'*R(Nuv:end)/(F0(Nuv:end)'*F0(Nuv:end))); rh=real(rh); 
        end

    end

    
    
end

