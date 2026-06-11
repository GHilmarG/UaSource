function ddIddCfp=Calc_FixPoint_ddIddC(CtrlVar,MUA,ub,vb,ud,vd,C,GF)
    
    % a rather add-hoc estimate of the Hessian of the data-misfit with respect to C
    % No garantee, but often dramatically improves rates of convergence.
    %
    
    
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
    % dI^2/dC^2   = d/du  ( -(umeas-u) (tb/He)^m)
    %             = -du/dC  (tb/He)^m                    (again ignoring the dependecey of tb on C)
    %             = -(tb/He)^(2m)
    %             = -(C/u)^2
    %
    %    on the floating part He->0 and tb->0 so what is the right expression for  the floating part?
    %    It appears plausible that dI/dC is then zero and ddIddC too, however the Hessian can not be zero
    %    as that would create infinity large Newton step for any non-zero righ-hand side
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
    
    if CtrlVar.CisElementBased
        ub=Nodes2EleMean(MUA.connectivity,ub); 
        vb=Nodes2EleMean(MUA.connectivity,vb); 
        
    end
    
    speed=sqrt(ub.*ub+vb.*vb);
    
    ddIddCfp=((speed+CtrlVar.SpeedZero)./(C+CtrlVar.Cmin)).^2;
    %ddIddCfp=GF.ele.*ddIddCfp;
    ddIddCfp=CtrlVar.MisfitMultiplier*ddIddCfp;
    
%     % make it less dependent on local conditions
%     Nnodes=length(u);
%     M= Ele2Nodes(connectivity,Nnodes) ;
%     ddIddCfp = M*ddIddCfp;
%     ddIddCfp=Nodes2EleMean(connectivity,ddIddCfp);
    
    % the calculated Hessian is always positive, by construction. I must make sure it does not approach zero 
    % so I define some plausible lower limit
    MinddIddCfp=median(ddIddCfp)/1e5; 
    ddIddCfp(ddIddCfp<MinddIddCfp)=MinddIddCfp;
    
    if CtrlVar.MisfitMultiplier==0
        ddIddCfp=ddIddCfp*0;
    end
    
    N=length(C);
    
    ddIddCfp=sparse(1:N,1:N,ddIddCfp);
    
end

