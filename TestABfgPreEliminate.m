
%%

A=[0 1 0 ; ...
   1 0 0 ;
   0 0 1 ];
B=[] ; 
f=[1 ; 1 ; 1 ];
g=[] ; 

CtrlVar=[]; 


[x,y]=ABfgPreEliminate(CtrlVar,A,B,f,g)

%%


A=[0 0 0 0; ...
   0 0 0 0; ...
   0 0 1 0; ...
   0 0 0 1 ];
B=[1 0  0 0 ; ...
   0 1  0 0  ] ;
f=[10 ; 2 ; 10 ; 1000];
g=[10  ; 10 ] ;

CtrlVar=[];

[nA,mA]=size(A) ; [nB,mB]=size(B) ; 
iW=zeros(nB,nB); 

fg=[f;g];

[x,y,tolA,tolB]=ABfgPreEliminate(CtrlVar,A,B,f,g);

T=[A B' ; B iW];
norm(T\fg-[x;y])

%%
CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=1;
CtrlVar.InfoLevelLinSolve=10;
CtrlVar.ALSIterationMin=1;
CtrlVar.ALSIterationMax=100;
CtrlVar.ALSpower=2;
CtrlVar.Solve.LUvector=true;
CtrlVar.TestForRealValues=true;
CtrlVar.LinSolveTol=1e-6;
[x,y] = AugmentedLagrangianSolver(A,B,f,g,[],CtrlVar);



%%
