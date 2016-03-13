
% both node definitions work
node=[-1 -1 ; 0 0  ; 1 -1 ; 0 1 ; -1 -1];
node=[-1 -1 ; 0 0  ; 1 -1 ; 0 1 ];

p=[0 0 ; -1 -1 ; 1 -1 ; 0 1] ;

figure
plot(node(:,1),node(:,2))


[in,on]=inpoly(p,node)

