function ddIregddC=Calc_ddIregddC(CtrlVar,MUA,CC)

narginchk(3,3)

%
% returns an approximation to the Hessian using only diagonal elements of CC
%
% The regularisation term is
%      Cres=(C-C_prior);
%      IRegC=CtrlVar.RegCMultiplier*Cres'*(CC\Cres)/(2N)    ;
%    dIregdC=CtrlVar.RegCMultiplier*(CC\(C-C_prior)) / N;
%  ddIregddC=CtrlVar.RegCMultiplier*inv(CC) /N

[N,M]=size(CC);

if CtrlVar.isRegC
    ddIregddC=CtrlVar.RegCMultiplier./spdiags(CC,0);
    ddIregddC=sparse(1:N,1:N,ddIregddC)/N;

    
else
    
    if CtrlVar.CisElementBased
        ddIregddC=sparse(MUA.Nele,MUA.Nele);
    else
        ddIregddC=sparse(MUA.Nnodes,MUA.Nnodes);
    end
end



end