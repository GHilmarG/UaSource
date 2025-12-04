function [Values,FA2B]=ExtrapolateFromNodesAtoNodesB(CtrlVar,x,y,NodesA,NodesB,Values)

% Extrapolates from the subset of nodes, NodesA, to another subset of nodes, NodesB, given the values ValuesA on the NodesA
%
% x & y are the nodal coordinates of the WHOLE mesh and not just the x, y values of NodesA.
%
% NodesA and NodesB should be logical arrays
%
% 'Extrapolation' can here involve either interpolation or extrapolation!
% 
% Values over the nodes in set B (NodesB) within the convex hull of nodesA are interpolated. Extrapolatin onto NodesB is only
% used for those nodes outside of the convex hull of NodesA.


if numel(x) ~= numel(Values)

    Values=[] ; FA2B=[] ;
    warning("ExtrapolateFromNodesAtoNodesB:IncorrectInputs","NodesA and ValuesA must have same number of elements!")
    return

end

if numel(find(NodesB))==0
    Values=[] ; FA2B=[] ;
    return
end


FA2B=scatteredInterpolant();
FA2B.Points=[x(NodesA) y(NodesA)];
FA2B.Values=Values(NodesA) ;

FA2B.ExtrapolationMethod='nearest'; FA2B.Method='linear';

Values(NodesB)=FA2B(x(NodesB),y(NodesB)) ;


end