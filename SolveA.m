function x=SolveA(A,b)
%#codgen

%% codegen -config:mex SolveA -args {zeros(10,10) , zeros(10,1)}
%% codegen -config:mex SolveA -args {coder.typeof(zeros(10,10),[inf inf]), coder.typeof(zeros(10,1),[inf 1])}
%% codegen -config:mex SolveA -args {coder.typeof(zeros(10,10),[10000 10000],1), coder.typeof(zeros(10,1),[10000 1],1)}
%% n=10000 ; A=rand(n,n) ; b=rand(n,1) ; tic ; x=SolveA_mex(A,b); toc


x=A\b;

end