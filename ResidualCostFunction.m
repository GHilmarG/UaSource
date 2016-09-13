function  [r,rl,ruv,rh]=ResidualCostFunction(R1,R2,F0,Nnodes)

nargoutchk(1,4)
narginchk(4,4)

% square of normalized residual
% returns (R1'*R1+R2'*R2)/(F0'*F0);
%
% ruv and rh are normalized with F0'F0
% rl is not normalized

if isempty(R2) ; R2=0 ; end

ruv=R1'*R1;
rl=R2'*R2;

r=full(real((ruv+rl)/(F0'*F0)));


if nargout > 2
    
    Nuv=2*Nnodes;
    
    ruv=full(R1(1:Nuv)'*R1(1:Nuv)/(F0(1:Nuv)'*F0(1:Nuv)));
    ruv=real(ruv);
    
    if nargout>3
        rh=full(R1(Nuv+1:end)'*R1(Nuv+1:end)/(F0(Nuv+1:end)'*F0(Nuv+1:end)));
        rh=real(rh);
    end
    
    
end



end

