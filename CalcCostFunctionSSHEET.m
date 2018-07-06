function r=CalcCostFunctionSSHEET(CtrlVar,gamma,dh,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt,Lh,lambdah,dlambdah,F0,ch)

s1=s1+gamma*dh;
h=s1-b1;
R=MatrixAssemblySSHEETtransient2HD(CtrlVar,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt);


if ~isempty(Lh)
    
    frhs=-R-Lh'*(lambdah+gamma*dlambdah);
    grhs=ch-Lh*(h+gamma*dh);

else
    frhs=-R;
    grhs=[];
end

r=ResidualCostFunction(frhs,grhs,F0,MUA.Nnodes);



end

