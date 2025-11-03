
function [ufixednode,ufixedvalue,vfixednode,vfixedvalue,...
    utiedA,utiedB,vtiedA,vtiedB,FixedNormalVelocityNode,FixedNormalVelocityValue,...
    hfixednode,hfixedvalue,htiedA,htiedB]=FromMLCtoBCs(Luv,Luvrhs,Lh,Lhrhs)

% Takes the multi-linear constraint matrices, on which the BCs are based internally
% and returns arrays of fixed nodal values, ties, and so on. 
% In some ways the inverse of `FromBCsMLC.m'
%
% Note: TBD needs updating



[ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,FixedNormalVelocityNode,FixedNormalVelocityValue]=...
    FromLuvTouvBCs(Luv,Luvrhs);

[hfixednode,hfixedvalue,htiedA,htiedB]=FromLhTohBCs(Lh,Lhrhs);

end


