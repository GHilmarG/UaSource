function Reactions=CalculateReactions(CtrlVar,MUA,BCs,l)

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
%   l is one of the outputs of Ua available in UaOutputs
%
%   If MLC is not available, calculate MLC using MLC=BCs2MLC(MUA,BCs)
%
%   Example:
%   To calculate and plot reactions from within UaOutputs
%   MLC=BCs2MLC(MUA,BCs) ;
%   Reactions=CalculateReactions(CtrlVar,MLC,l)
%   PlotReactions(CtrlVar,MUA,Reactions);
%
% Reactions are defined for all the nodes, but for nodes where no BCs have been applied,
% they will automatically be equal to zero. However, in the special case where no essential
% BCs are applied, reactions are returned as an empty matrix.
%
%
   
MLC=BCs2MLC(CtrlVar,MUA,BCs) ;

if ~CtrlVar.LinFEbasis
    if ~isfield(MUA,'M')
        MUA.M=MassMatrix2D1dof(MUA);
    end
    
    Mblock=MassMatrixBlockDiagonal2D(MUA);
end

if ~isempty(l.ubvb)
    if CtrlVar.LinFEbasis
        Reactions.ubvb=MLC.ubvbL'*l.ubvb;
    else
        luv=MLC.ubvbL'*l.ubvb;
        Reactions.ubvb=MUA.M\(luv(1:MUA.Nnodes)+luv(MUA.Nnodes+1:end)); 
    end
else
    Reactions.ubvb=[];
end

if ~isempty(l.udvd)
    if CtrlVar.LinFEbasis
        Reactions.udvd=MLC.udvdL'*l.udvd;
    else
        Reactions.udvd=Mblock\MLC.udvdL'*l.udvd;
    end
else
    Reactions.udvd=[];
end

if ~isempty(l.h)
    if CtrlVar.LinFEbasis
        Reactions.h=MLC.hL'*l.h;
    else
        Reactions.h=MUA.M\MLC.hL'*l.h;
    end
else
    Reactions.h=[];
end


end

