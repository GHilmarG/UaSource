function Mblock=MassMatrixBlockDiagonal2D(MUA)


% calculates the mass matrix, ie : int N_p N_q,

if ~isfield(MUA,'M')
    MUA.M=MassMatrix2D1dof(MUA);
end


Z=sparse(MUA.Nnodes,MUA.Nnodes);
Mblock=[ MUA.M Z ; Z MUA.M ];


end



