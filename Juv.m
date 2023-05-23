
function [Juv,dJduv,Huv,RunInfo]=Juv(uv,UserVar,RunInfo,CtrlVar,MUA,F,fext0,Aeq,beq)



N=MUA.Nnodes;

F.ub=uv(1:N);
F.vb=uv(N+1:2*N);


if nargout>2

    [Ruv,Kuv]=KRTFgeneralBCs(CtrlVar,MUA,F);

else

    Ruv=KRTFgeneralBCs(CtrlVar,MUA,F);
    Kuv=[];

end

if ~isempty(L)
    f=-Ruv-L'*l;
    g=beq-Aeq*[F.ub;F.vb];

else
    f=-Ruv;
    g=[];
    % dl=[];
end



% Juv=full((f'*f)./(fext0'*fext0+1000*eps));

Juv=full([f;g]'*[f;g]./(fext0'*fext0+1000*eps));

% Mblock=MassMatrixBlockDiagonal2D(MUA); M=Mblock;


dJduv= Ruv;



if nargout>2
    Huv=Kuv;
    fprintf("Hessian \n")
end

%%



end



