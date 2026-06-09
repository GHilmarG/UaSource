
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

if CtrlVar.BCsRowSubsetSelection

    fprintf("AssembleLuvSSTREAM: Checking if uv boundary conditions are linearly independent.\n")
    [Luv,row_idx,flag] = RowSubsetSelection(Luv);

    if flag==1
        cuv=cuv(row_idx);
        fprintf("AssembleLuvSSTREAM: uv boundary conditions were found NOT to be linearly independent.\n")
        fprintf("This has now been corrected internally using row-subset selection, but it might be good to reconsider how the BCs were defined in DefineBoundaryConditions.m \n")
    end

end

end

