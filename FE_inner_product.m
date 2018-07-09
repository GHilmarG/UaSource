function [Prod,M]=FE_inner_product(f,g,M)

% [Prod,M]=FE_inner_product(f,g,MUA,M)
%  calulates <f|M|g> where M is the mass matrix
% M=MassMatrix2D1dof(MUA);


Prod=f'*M*g;

end



