
%%

load solveKApePIGTWGuvh250896.mat

%%

CtrlVar.AsymmSolver="auto";


tic
[x,y]=solveKApe(A,B,f,g,x0,y0,CtrlVar) ; 
toc


%%

CtrlVar.InfoLevelLinSolve=100; 
CtrlVar.AsymmSolver='EliminateBCsSolveSystemDirectly';


CtrlVar.AsymmSolver='EliminateBCsSolveSystemIterativly';

tic
[x,y]=solveKApe(A,B,f,g,x0,y0,CtrlVar) ; 
toc
