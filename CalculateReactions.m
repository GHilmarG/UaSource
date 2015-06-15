function Reactions=CalculateReactions(MLC,l)

% Reactions=CalculateReactions(MLC,l)
% calculates nodal reactions
% MLC  : muliple-linear-constraint matrix 
%   l  : Lagrange variables 
%
%   l is one of the outputs of Ua available in UaOutputs
%
%   If MLC is not availabel, calculate MLC using MLC=BCs2MLC(MUA,BCs)
%
%   Example: 
%   To calculate and plot reactions from within UaOutputs
%   MLC=BCs2MLC(MUA,BCs) ; 
%   Reactions=CalculateReactions(MLC,l)
%   PlotReactions(CtrlVar,MUA,Reactions);
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

