function DriveSolveA

n=100;

A=rand(n,n); b=rand(n,1);

tic
x=SolveA(A,b);
toc

%tic
%x=SolveA_mex(A,b);
%toc

end

