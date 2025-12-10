load PIG-TWG-RestartFile.mat

CtrlVarInRestartFile.MUA.MassMatrix=true;
MUA=UpdateMUA(CtrlVarInRestartFile,MUA);


x1=F.s ; x2=F.s ;
tic ;

N=100; 
for ii=1:N
    x1=MUA.M\x1 ;
    x1=x1/norm(x1);
end
tM=toc;

tic
dM=decomposition(MUA.M,"chol") ; 

for ii=1:N
    x2=dM\x2 ;
    x2=x2/norm(x2);
end
tdM=toc;

fprintf("\n\n    time with  \t \t without   \t\t speed up \n")
fprintf("\t %f \t\t %f \t\t %f\n",tM,tdM,tM/tdM)