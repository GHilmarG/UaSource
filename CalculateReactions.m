function Reactions=CalculateReactions(CtrlVar,MUA,BCs,l)

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
     if ~isfield(MUA,'dM')
        MUA.dM=decomposition(MUA.M);
    end
   
end

if ~isempty(l.ubvb)
    if CtrlVar.LinFEbasis
        Reactions.ubvb=MLC.ubvbL'*l.ubvb;
    else
        luv=MLC.ubvbL'*l.ubvb;
        Rx=MUA.dM\luv(1:MUA.Nnodes);
        Ry=MUA.dM\luv(MUA.Nnodes+1:end); 
        Reactions.ubvb=[Rx;Ry];
    end
else
    Reactions.ubvb=[];
end

if ~isempty(l.udvd)
    if CtrlVar.LinFEbasis
        Reactions.udvd=MLC.udvdL'*l.udvd;
    else
        luv=MLC.udvdL'*l.udvd;
        Reactions.udvd(1:MUA.Nnodes)=MUA.dM\luv(1:MUA.Nnodes);
        Reactions.udvd(MUA.Nnodes+1:end)=MUA.dM\luv(MUA.Nnodes+1:end);

    end
else
    Reactions.udvd=[];
end

if ~isempty(l.h)
    if CtrlVar.LinFEbasis
        Reactions.h=MLC.hL'*l.h;
    else
        Reactions.h=MUA.dM\(MLC.hL'*l.h);
    end
else
    Reactions.h=[];
end


end

