%%
coordinates=[ 1 0 ; 0 0 ; 0 1]; connectivity=[ 1 2 3];

fxx=[1 ; 1 ; 1 ] ;  fyy=[1 ; 1; 1]; fxy=[0 ; 0 ; 0 ] ;  fyx=[0 ; 0; 0];
bx=[1 ; 1 ; 1 ] ;  by=[1 ; 1; 1];


%%
% 
 coordinates=[ 1 0 ; 0 0 ; 0 1 ; 1 1 ; 0.5 0.5 ]; connectivity=[ 1 2 5 ; 2 3 5 ; 3 4 5 ; 4 1 5];
 
 fxx=[1 ; 1 ; 1 ; 1 ; 1] ;  fyy=[1 ; 1; 1 ; 1 ; 1];
 fxy=[0 ; 0 ; 0 ; 0 ; 0] ;  fyx=[0 ; 0; 0 ; 0 ; 0];
 bx=[1 ; 1 ; 1 ; 1 ; 1] ;  by=[1 ; 1; 1 ; 1 ; 1];

%%

CtrlVar.TriNodes=3;

[K,rhs]=FE_BoundaryIntegral(CtrlVar,coordinates,connectivity,fxx,fyy,fxy,fyx,bx,by);


