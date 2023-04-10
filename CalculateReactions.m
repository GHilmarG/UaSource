function [Reactions,lStar]=CalculateReactions(CtrlVar,MUA,BCs,l)

% save TestSaveCalculateReactions
narginchk(4,4)

%%
%
%   Reactions=CalculateReactions(MLC,l)
%
% calculates nodal reactions
% Reactions=CalculateReactions(MLC,l)
%
% MLC  : muliple-linear-constraint matrix
%   l  : Lagrange variables
%
%   l is one of the outputs of Ua available in DefineOutputs
%
%   If MLC is not available, calculate MLC using MLC=BCs2MLC(MUA,BCs)
%
%   Example:
%   To calculate and plot reactions from within DefineOutputs
%   MLC=BCs2MLC(MUA,BCs) ;
%   Reactions=CalculateReactions(CtrlVar,MLC,l)
%   PlotReactions(CtrlVar,MUA,Reactions);
%
% Reactions are defined for all the nodes, but for nodes where no BCs have been applied,
% they will automatically be equal to zero. However, in the special case where no essential
% BCs are applied, reactions are returned as an empty matrix.
%
% If the multi-linear contraint matrix L was assembled in the FE basis, then the reactions
% are already physically correct and I only need to split them up in uv and h reactions.
%
% If L is a point-contraint matrix (default) then I need to map these into physical space
% usuing
%
%   lambda^* = M^{-1} L' lambda 
%
%  This gives lamnda over all nodes.
%
% To restrict it to the nodes over which the constrants were applied use:
%
%     lambda^* = (L L')^{-1} L M^{-1} L' lambda 
%
%
%% 

Reactions.ubvb=[];
Reactions.udvd=[];
Reactions.h=[];
lStar.h=[];  % only calculating physical lambdas for thickness constraints

MLC=BCs2MLC(CtrlVar,MUA,BCs) ;

if ~CtrlVar.LinFEbasis
    if ~isfield(MUA,'M')
        MUA.M=MassMatrix2D1dof(MUA);
    end
end

if ~isempty(l.ubvb)
    
        luv=MLC.ubvbL'*l.ubvb;
        Rx=MUA.M\luv(1:MUA.Nnodes);
        Ry=MUA.M\luv(MUA.Nnodes+1:end); 
        Reactions.ubvb=full([Rx;Ry]);
        
        
end

if ~isempty(l.udvd)
    
    luv=MLC.udvdL'*l.udvd;
    Reactions.udvd(1:MUA.Nnodes)=full(MUA.M\luv(1:MUA.Nnodes));
    Reactions.udvd(MUA.Nnodes+1:end)=full(MUA.M\luv(MUA.Nnodes+1:end));
end

if ~isempty(l.h)
    M=MUA.M;
    L=MLC.hL;
    lambda=l.h ;
    Reactions.h=full(M\(L'*lambda));
    
    if nargout>1
        lStar.h=full((L*L')\(L*Reactions.h));
    end
    
end




end

