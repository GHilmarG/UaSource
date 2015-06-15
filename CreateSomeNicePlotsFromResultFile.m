load PIGDumpFile100

%%
GF = GL2d(B,S,h,rhow,rho,connectivity,CtrlVar);
H0=1000;
Resolution=500;
x=coordinates(:,1); y=coordinates(:,2);

xgrid= min(DTxy.X(:,1)):Resolution:max(DTxy.X(:,1)) ;
ygrid= min(DTxy.X(:,2)):Resolution:max(DTxy.X(:,2)) ;
[X,Y]=meshgrid(xgrid,ygrid);


%d=Draft(S,B,b,CtrlVar);
%[Draftx,Drafty] = FindGL(DTxy,d,Resolution);

BoundaryNodes=FindBoundaryNodes(connectivity,coordinates);
GF.node(BoundaryNodes)=0;
[GLx,GLy,GLxUpper,GLyUpper,GLxLower,GLyLower] = FindGL(DTxy,GF.node,Resolution);
    
%%


xp=100 ; yp=100;figure('Position',[xp yp xp+400 yp+400]) ; 
trisurf(TRIxy,x/H0,y/H0,log10(C),'EdgeColor','none') ;  title('log10(Cest)') ;lightangle(-45,30) ; lighting phong ; 
Cbar=colorbar ; title(Cbar,'(m a^{-1} kPa^{-m})')
view(0,90) ; 
axis tight ; xlabel('xps (km)') ; ylabel('yps (km)') ; title('log_{10}(Basal slipperiness)')
aspect=get(gca,'DataAspectRatio') ; set(gca,'DataAspectRatio',[aspect(1) aspect(1) aspect(3)]) ; 
%axis([-1700 -1500 -410 -100]) ; 
hold on; plot3(GLx/H0,GLy/H0,GLx*0+max(C),'r','LineWidth',1)


xp=xp+50 ; yp=yp+50;figure('Position',[xp yp xp+400 yp+400]) ; 
trisurf(TRIxy,x/H0,y/H0,AGlen,'EdgeColor','none') ;  title('log10(Cest)') ;lightangle(-45,30) ; lighting phong ; 
Cbar=colorbar ; title(Cbar,'(a^{-1} kPa^{-3})')
hold on; plot3(GLx/H0,GLy/H0,GLx*0+max(AGlen),'r','LineWidth',1)
view(0,90) ; 
axis tight ; xlabel('xps (km)') ; ylabel('yps (km)') ; title('Rate factor') ; set(gca,'DataAspectRatio',[1e10 1d10 1]) ; axis([-1700 -1500 -410 -100]) ; 

%%

xp=50 ; yp=50;figure('Position',[xp yp xp+500 yp+500]) ;  
speed=sqrt(u.*u+v.*v);
trisurf(TRIxy,x/H0,y/H0,speed,'EdgeColor','none') ;  title('speed') ;lightangle(-45,30) ; lighting phong ; 
Cbar=colorbar ; title(Cbar,'(m a^{-1})')
hold on; plot3(GLx/H0,GLy/H0,GLx*0+max(speed),'r','LineWidth',1)
view(0,90) ; 
aspect=get(gca,'DataAspectRatio') ; set(gca,'DataAspectRatio',[aspect(1) aspect(1) aspect(3)]) ; 
axis tight ; xlabel('xps (km)') ; ylabel('yps (km)') ; title('Speed') ; %set(gca,'DataAspectRatio',[1e10 1d10 1]) ; axis([-1700 -1500 -410 -100]) ; 

%%
speed=sqrt(uMeas.*uMeas+vMeas.*vMeas);
xp=50 ; yp=50;figure('Position',[xp yp xp+500 yp+500]) ;  
trisurf(TRIxy,x/H0,y/H0,speed,'EdgeColor','none') ;  title('speed') ;lightangle(-45,30) ; lighting phong ; 
Cbar=colorbar ; title(Cbar,'(m a^{-1})')
hold on; plot3(GLx/H0,GLy/H0,GLx*0+max(speed),'r','LineWidth',1)
view(0,90) ; 
aspect=get(gca,'DataAspectRatio') ; set(gca,'DataAspectRatio',[aspect(1) aspect(1) aspect(3)]) ; 
axis tight ; xlabel('xps (km)') ; ylabel('yps (km)') ; title('Measured Speed') ; %set(gca,'DataAspectRatio',[1e10 1d10 1]) ; axis([-1700 -1500 -410 -100]) ; 
%%
