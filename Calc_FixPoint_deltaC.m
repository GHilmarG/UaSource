function dIdC=Calc_FixPoint_deltaC(CtrlVar,MUA,C,m,GF,ub,vb,usMeas,vsMeas)
    
    % Fix point idea:
    %
    % u=C tb^m   ; minimize (umeas-u) = umeas-C tb^m by setting umeas-u=0 using Cnew=umeas/tb^m
    % if I write Cnew=C + delta C = C+(Cnew-C) then delta C= C umeas/u - C = (umeas/u-1)C
    % or delta C= C (umeas-u)/u
    %
    % tb= He (u/c)^(1/m) -> u= C (tb/He)^m
    % ignore the dependency of tb on C
    % u_meas=Ctrue (tb/He)^m
    % u_num= C (tb/He)^m
    % Ctrue=C+(Ctrue-C)=C+deltaC=C+(u_meas-u_num) (tb/He)^{-m}
    % I want u=u_meas  ->  deltaC=(u_meas-u0) (tb/He)^{-m}
    %                            =(u_meas-u0) C/speed
    %
    % Compare to NR applied to cost function:
    % I= 1/2 (umeas-u)^2
    % dI/dC=dI/du du/dC
    %       =(umeas-u) du/dC
    %       =(umeas-u) (tb/He)^m   (if I assume that dtb/dC=0 then du/dC=(tb/He)^m)
    %
    % dI^2/dC^2 C = -du/dC  (tb/He)^m
    %             = -(tb/He)^(2m)
    %
    % The Newton system is:  dI^2/dC^2 Delta C = - dI/dC
    %                      -(tb/He)^(2m) Delta C= (umeas-u) (tb/He)^m
    %                                    Delta C= -(umeas-u) (He/tb)^{m}
    %
    %   Delta C= -H\dI/dC= (He/tb)^(m)  (same result as from fix-point idea)
    %           =-(umeas-u) C/u
    %
    %
    %  tb=He .* C.^(-1/m).*speed.^(1/m);
    %  (tb/He)^m=  speed/C and (He/tb)^m=C/speed
    %
    % One problem with using Delta C = (umeas-u) (tb/He)^(-m)
    %                                = (umeas-u) C/speed
    % is that tb u is not dependent on C where floating
    % we can think about this as if C is then infinitly large
    % and that the actual slipperiness is C/He
    % Hence: Delta C= He (umeas-u) C/speed
    %
    
    
    if CtrlVar.CisElementBased
        uEle=Nodes2EleMean(MUA.connectivity,ub); 
        vEle=Nodes2EleMean(MUA.connectivity,vb); 
        speed=sqrt(uEle.*uEle+vEle.*vEle);
        uEleMeas=Nodes2EleMean(MUA.connectivity,usMeas); vEleMeas=Nodes2EleMean(MUA.connectivity,vsMeas);
        speedMeas=sqrt(uEleMeas.*uEleMeas+vEleMeas.*vEleMeas);
    else 
        speed=sqrt(ub.*ub+vb.*vb);
        speedMeas=sqrt(usMeas.*usMeas+vsMeas.*vsMeas);
    end
    
    
    %dIdC=-(speedEleMeas-speedEle)./(tb./GF.ele).^m; % this is not based on dIdC= dI/du  du/dC, this is an estimate of the (negative) Newton step -H\grad I
    
    minSpeed=10;
    speedMeas(speedMeas<minSpeed)=minSpeed; % see Remark
    dIdC=-C.*(speedMeas-speed)./(speed+minSpeed) ; % this is not based on dIdC= dI/du  du/dC,
    % but an estimate of the (negative) Newton step -H\grad I
    % the problem with this expression is that u is not equal to c (tb/He)^m when He->0
    
    % Remark: when speedEleMeas=0 then  dIdC=-(speedEleMeas-speedEle)./(speedEle./C)=1/C , irrespectivly of how close speedEle is to zero.
    % This will drive C further and further towards zero with every iteration.
    % in some sense this is correct, because only for C stricly equal to 0 is speedEle=0
    % The problme is that this creates infinitly large beta^2
    % and to avoid this I introduce a minimum speed of 1 m/a
    
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        dIdC=log(10)*C.*dIdC;
    end
    
    
    if CtrlVar.CisElementBased
        dIdC=dIdC.*GF.ele;
    else
        dIdC=dIdC.*GF.node;
    end
    
   
    
    
    
    
end

