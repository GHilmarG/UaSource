


load TestSave DTxy DTint MeshBoundaryNodes coordinates b B s h u v w C

x=coordinates(:,1) ; y=coordinates(:,2);

TRIxy=DTxy.Triangulation;
TRIint=DTint.Triangulation;

%%
ic=incenters(DTxy);

[cn,on] = inpoly(ic,MeshBoundaryNodes);

%%


figure ; trisurf(TRIxy(cn,:),x,y,b) ;  title('b') ;xlabel('x') ; ylabel('y')

%%
