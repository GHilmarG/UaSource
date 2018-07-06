function M=MassMatrixBlockDiagonal2D(MUA)


% calculates the mass matrix, ie : int N_p N_q,

M=MassMatrix2D1dof(MUA);

Z=sparse(MUA.Nnodes,MUA.Nnodes);
M=[ M Z ; Z M ];


end



