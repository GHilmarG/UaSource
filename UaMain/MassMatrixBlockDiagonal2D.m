function Mblock=MassMatrixBlockDiagonal2D(MUA)


% calculates the mass matrix, ie : int N_p N_q,

if ~isfield(MUA,'M')
    MUA.M=MassMatrix2D1dof(MUA);
end

% Note: To do : pretty sure I could use blkdiag(MUA.M,MUA.m) here instead.

Z=sparse(MUA.Nnodes,MUA.Nnodes);
Mblock=[ MUA.M Z ; Z MUA.M ];


end



