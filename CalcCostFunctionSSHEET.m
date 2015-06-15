function r=CalcCostFunctionSSHEET(CtrlVar,gamma,dh,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt,Lh,lambdah,dlambdah,F0)

s1=s1+gamma*dh;
R=MatrixAssemblySSHEETtransient2HD(CtrlVar,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt);

if ~isempty(lambdah)
    R=R+Lh'*(lambdah+gamma*dlambdah);
end

r=ResidualCostFunction(R,F0);



end

