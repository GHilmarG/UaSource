function Reactions=CalculateReactions(MLC,l)

% calculates nodal reactions
% Reactions=CalculateReactions(MLC,l)
% 
% MLC  : muliple-linear-constraint matrix 
%   l  : Lagrange variables 
%
%   l is one of the outputs of Ua available in UaOutputs
%
%   If MLC is not available, calculate MLC using MLC=BCs2MLC(MUA,BCs)
%
%   Example: 
%   To calculate and plot reactions from within UaOutputs
%   MLC=BCs2MLC(MUA,BCs) ; 
%   Reactions=CalculateReactions(MLC,l)
%   PlotReactions(CtrlVar,MUA,Reactions);
%
% Reactions are defined for all the nodes, but for nodes where no BCs have been applied, 
% they will automatically be equal to zero. However, in the special case where no (non-natural) 
% BCs are applied, reactions are returned as an empty matrix.
% 
%

if ~isempty(l.ubvb)
    Reactions.ubvb=MLC.ubvbL'*l.ubvb;
else
    Reactions.ubvb=[];
end

if ~isempty(l.udvd)
    Reactions.udvd=MLC.udvdL'*l.udvd;
else
    Reactions.udvd=[];
end

if ~isempty(l.h)
    Reactions.h=MLC.hL'*l.h;
else
    Reactions.h=[];
end


end

