function ddIregddC=Calc_ddIregddC(CtrlVar,CC)
    
    %
    % returns an approximation to the Hessian using only diagonal elements of CC
    %
    % The regularisation term is 
    %      Cres=(C-C_prior);
    %      IRegC=CtrlVar.RegCMultiplier*Cres'*(CC\Cres)/2    ;
    %    dIregdC=CtrlVar.RegCMultiplier*(CC\(C-C_prior));
    %  ddIregddC=CtrlVar.RegCMultiplier*inv(CC)
    
    [N,M]=size(CC);
    
     if CtrlVar.isRegC
        ddIregddC=CtrlVar.RegCMultiplier./spdiags(CC,0);
    else
        ddIregddC=zeros(N,1);
     end
     
     ddIregddC=sparse(1:N,1:N,ddIregddC);
        
end