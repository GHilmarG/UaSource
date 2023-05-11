
function [Juv,dJduv,Huv,RunInfo]=Juv(uv,UserVar,RunInfo,CtrlVar,MUA,F,fext0)



N=MUA.Nnodes;

F.ub=uv(1:N);
F.vb=uv(N+1:2*N);


if nargout>2

    [Ruv,Kuv]=KRTFgeneralBCs(CtrlVar,MUA,F);

else

    Ruv=KRTFgeneralBCs(CtrlVar,MUA,F);
    Kuv=[];

end



f=-Ruv;
Juv=full((f'*f)./(fext0'*fext0+1000*eps));
dJduv=Ruv;
Huv=Kuv;


end



