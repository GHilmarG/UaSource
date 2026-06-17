
function [Luv,cuv]=AssembleLuvSSTREAM(CtrlVar,MUA,BCs)

MLC=BCs2MLC(CtrlVar,MUA,BCs);
Luv=MLC.ubvbL;
cuv=MLC.ubvbRhs;


% Luv  : #uv constrains x 2MUA.Nnodes
% L=[Luv]



if CtrlVar.LinFEbasis && ( numel(Luv)>0 || numel(cuv)>0 )
    
    % L M L' L
    %
    % L -> L M
    % c -> L M L' c
    if ~isfield(MUA,'M')
        MUA.M=MassMatrix2D1dof(MUA);
    end
    
    Mblock=MassMatrixBlockDiagonal2D(MUA);
    
    if numel(Luv)>0
        Luv=Luv*Mblock ;
    end
 
    
    if numel(cuv)>0
        cuv=(Luv*Mblock*Luv')*cuv ;
    end
    
 
end


end

