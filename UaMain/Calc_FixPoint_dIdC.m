function dIdCfp=Calc_FixPoint_dIdC(CtrlVar,MUA,C,m,GF,u,v,uMeas,vMeas,Cd)
    
 
    % ignore the dependency of tb on C

    % u= C (tb/He)^m
    %
    % Compare to NR applied to cost function:
    % I= 1/2 (umeas-u)^2
    % dI/dC=dI/du du/dC
    %       =-(umeas-u) du/dC
    %       =-(umeas-u) (tb/He)^m   (if I assume that dtb/dC=0 then du/dC=(tb/He)^m)
    %       =-(umeas-u) u/C
    %
    % dI^2/dC^2 C = du/dC  (tb/He)^m
    %             = (tb/He)^(2m)
    %             = (C/u)^2
    %
    % The Newton system is:  dI^2/dC^2 Delta C = - dI/dC
    %                      -(tb/He)^(2m) Delta C= (umeas-u) (tb/He)^m
    % or                        -(u/C)^2 Delta C= (umeas-u) u/C
    %                                    Delta C= -(umeas-u) (He/tb)^{m}
    %                                           = -(umeas-u) C/u
    %   Delta C= -H\dI/dC= (He/tb)^(m)  (same result as from fix-point idea)
    %           =-(umeas-u) C/u
    %
    % Hence: Delta C= He (umeas-u) speed/C
    %
    %Cd=sparse(1:2*Nnodes,1:2*Nnodes,[uError.^2;vError.^2],2*Nnodes,2*Nnodes);   % I=(B u-d)' inv(M) (B u -d )
    
    Nnodes=numel(u);
    DataErrors=diag(Cd,0);
    speedError2=DataErrors(1:Nnodes)+DataErrors(Nnodes+1:end);
    speedErrorEle2=Nodes2EleMean(MUA.connectivity,speedError2); 
    
    uEle=Nodes2EleMean(MUA.connectivity,u); vEle=Nodes2EleMean(MUA.connectivity,v); speedEle=sqrt(uEle.*uEle+vEle.*vEle);
    uEleMeas=Nodes2EleMean(MUA.connectivity,uMeas); vEleMeas=Nodes2EleMean(MUA.connectivity,vMeas);
    speedEleMeas=sqrt(uEleMeas.*uEleMeas+vEleMeas.*vEleMeas);
    minSpeed=10;
    
    dIdCfp=-(speedEleMeas-speedEle).*(speedEle+minSpeed)./C./speedErrorEle2;
    
    dIdCfp=dIdCfp.*GF.ele;
    
    dIdCfp=CtrlVar.MisfitMultiplier*dIdCfp;
    
    
end

