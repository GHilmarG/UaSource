function [Reactions,lStar]=CalculateReactions(CtrlVar,MUA,BCs,l,options)



%%
%

% Calculates nodal reactions.
%
%   Reactions=CalculateReactions(MLC,l)
%
% Nodal reactions are due to boundary conditions. They are here calculated from the values of the Lagrange parameters (l) introduced
% when enforcing Dirichlet boundary conditions.
%
% Reactions is a structure with the fields:
%
%   Reactions.ubvb
%   Reactions.udvd
%   Reactions.h
%
% By setting the optional argument 
%
%    hReactionsExcludeActiveSet=true  
%
%
% any positive thickness constraints that were added automatically, ie not by the user, during the run are excluded. 
%
% Note that, for example, the h reactions are related to the additional mass per area required to keep the thickness at
% prescribed values. This fictitious mass flux (a_h) can be calculated as
%
%   ah=Reactions.h/(F.rho*F.time)
%
% So if, for example, rho and time are in IS units, the units of the reactions are kg/m^2, and that of ah is m/yr.
%
% *Sign Convention:*  Same sign convention is used as when introducing the Lagrange parameters. It follows that if, for
% example, ADDITONAL mass flux is needed to keep the thickness at a given value, the resulting thickness reactions
% (Reactions.h) are NEGATIVE. This sign convention makes sense if active constraints are identified by the sign of the
% corresponding Lagrange parameters being negative. But this sign is presumably the opposite of what one might otherwise
% expect. So if the thickness reactions are NEGATIVE, then one needs to add a POSITIVE surface mass balance equal to
% 
%   -Reactions.h/(F.rho*F.time) 
% 
% to keep the thickness at desired value.
%
%
%
% l  : Lagrange variables
%
% l is one of the outputs of Ua available in DefineOutputs
%
%
%   Example:
%
% To calculate and plot reactions from within DefineOutputs
%
%   Reactions=CalculateReactions(CtrlVar,MUA,BCs,l); 
%   PlotReactions(CtrlVar,MUA,F,Reactions);
%
% Reactions are defined for all the nodes, but for nodes where no BCs have been applied,
% they will automatically be equal to zero. However, in the special case where no essential
% BCs are applied, reactions are returned as an empty matrix.
%
% If the multi-linear constraint matrix L was assembled in the FE basis, then the reactions
% are already physically correct and I only need to split them up in uv and h reactions.
%
% If L is a point-constraint matrix (default) then I need to map these into physical space
% using
%
%   lambda^* = M^{-1} L' lambda 
%
%  This gives lambda over all nodes.
%
% To restrict it to the nodes over which the constraints were applied use:
%
%     lambda^* = (L L')^{-1} L M^{-1} L' lambda 
%
%
%
% Note: If the active set has only been updated, but never applied, the dimensions of the Lagrange vector  (l.h) will be inconsistent with
%       the number of positive thickness constraints (BCs.hFixeNode).
%
% Also:  To plot reactions use:
%
%   PlotReactions(CtrlVar,MUA,F,Reactions)
%
%% 


arguments

    CtrlVar struct
    MUA     struct
    BCs     {mustBeA(BCs,{'struct','BoundaryConditions'})}
    l       {mustBeA(l,{'struct','UaLagrangeVariables'})}
    options.hReactionsExcludeActiveSet  logical = false   % this allows for h reactions to exclude those of the active thickness constraints
                                                            % only the reactions due to the user-defined thickness constraints are then calculated

end

if  ~isfield(MUA,"M") || isempty(MUA.M)
    MUA.M=MassMatrix2D1dof(MUA);
    MUA.dM=decomposition(MUA.M,'chol','upper') ;
end

if  ~isfield(MUA,"dM") || isempty(MUA.dM)
    MUA.dM=decomposition(MUA.M,'chol','upper') ;
end


Reactions.ubvb=[];
Reactions.udvd=[];
Reactions.h=[];
lStar.h=[];  % only calculating physical lambdas for thickness constraints

% getting rid of all positive thickness constraints as part of the active set
if options.hReactionsExcludeActiveSet
    BCs.hPosNode=[]; BCs.hPosValue=[]; lh=l.h(1:numel(BCs.hFixedNode)) ; l.h=lh ;
end


MLC=BCs2MLC(CtrlVar,MUA,BCs) ;

if ~CtrlVar.LinFEbasis
    if ~isfield(MUA,'M')
        MUA.M=MassMatrix2D1dof(MUA);
    end
end

if MUA.dM.MatrixSize(1) == 0   || MUA.dM.MatrixSize(2) == 0

    MUA=UpdateMUA(CtrlVar,MUA);

end

if ~isempty(l.ubvb)

    luv=MLC.ubvbL'*l.ubvb;

    try
        Rx=MUA.dM\luv(1:MUA.Nnodes);  % I can't see how I can know if the decomposition is in a valid state without just trying
    catch
        MUA.dM=decomposition(MUA.M,'chol','upper') ;
        Rx=MUA.dM\luv(1:MUA.Nnodes);
    end


    Ry=MUA.dM\luv(MUA.Nnodes+1:end);
    Reactions.ubvb=full([Rx;Ry]);


end

if ~isempty(l.udvd)
    
    luv=MLC.udvdL'*l.udvd;
    Reactions.udvd(1:MUA.Nnodes)=full(MUA.M\luv(1:MUA.Nnodes));
    Reactions.udvd(MUA.Nnodes+1:end)=full(MUA.M\luv(MUA.Nnodes+1:end));
end

if ~isempty(l.h)
    
    Lh=MLC.hL;
  
    Reactions.h=full(MUA.dM\(Lh'*l.h));
    
    if nargout>1
        lStar.h=full((Lh*Lh')\(Lh*Reactions.h));
    end
    
end




end

