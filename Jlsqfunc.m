




function  J=Jlsqfunc(CtrlVar,gamma,dx,dlambda,fun,L,c,x0,lambda0)



isLSQ=CtrlVar.lsqUa.isLSQ ;
CostMeasure=CtrlVar.lsqUa.CostMeasure;

x=x0+gamma*dx;
lambda0=lambda0+gamma*dlambda ;


if ~isempty(L)
    LTlambda=L'*lambda0 ;
    h =- (L*x-c);
else
    LTlambda=0;
    h=[];
end


if isLSQ
    [R,K]=fun(x) ;
    g =- (2*K'*R + LTlambda) ;
else
    R=fun(x) ;
    g =- (R + LTlambda) ;
end



if CostMeasure=="R2"

    J=full(R'*R) ;

elseif CostMeasure=="r2"

    d=[g;h]; 
    J=full(d'*d);


end

end




