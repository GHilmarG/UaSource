
function FE2dPlots(CtrlVar,DTxy,TRIxy,MeshBoundaryCoordinates,GF,dGFdt,coordinates,connectivity,...
        b,B,S,s,h,u,v,wSurf,dhdt,C,AGlen,m,n,xint,yint,wSurfInt,etaInt,exx,eyy,exy,e,time,...
        rho,rhow,a,as,ab,txx,tyy,txy)
    
fprintf(' FE2dPlots time=%-f \n ',time)

set(gcf, 'PaperType','A3') ; orient landscape ;  % usually better default values

%%
if ~isfield(CtrlVar,'PlotXYscale')
    CtrlVar.PlotXYscale=1000;
else
    CtrlVar.PlotXYscale=1000;  % take this out later
end


CtrlVar.GLresolutionWhenPlotting=1000; 
CtrlVar.GLtension=1e-11;

[Boundary.Nodes,Boundary.EdgeCornerNodes]=FindBoundaryNodes(connectivity,coordinates);  
xBoundary=coordinates(Boundary.EdgeCornerNodes,1);  yBoundary=coordinates(Boundary.EdgeCornerNodes,2); 

if ~all(GF.node==1)
    
    [GLx,GLy,GLxUpper,GLyUpper,GLxLower,GLyLower] = FindGL(DTxy,GF.node,CtrlVar,xBoundary,yBoundary);
    
    % For defining the coast  I need to `close' the grounding line, so I set all boundary nodes to floating status
    GFtemp=GF; GFtemp.node(Boundary.Nodes)=0;  
    [GLgeo,GLinfo]=GLgeometry(connectivity,coordinates,GFtemp,CtrlVar);
    [xGL,yGL,GLindex] = ArrangeGroundingLinePos(CtrlVar,GLgeo,1) ; % this is my definition of the `coast'
        
    % no I calculate all grounding lines element-wise
    [GLgeo,GLinfo]=GLgeometry(connectivity,coordinates,GF,CtrlVar);
    
    [xglc,yglc,nxGL,nyGL]=FindNiceLookingGLforPlottingPurposes('Contouring',DTxy,GF,CtrlVar,xBoundary,yBoundary);

    
    figure
    plot(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,'b') ; hold on
    plot(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,'g') ;
    plot(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,'r') ;
    plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
    
    axis equal
    title('Three different ways of calculationg GL pos')
    legend('Contouring of GF mask','Contouring GF mask followed by spline','coast','Elementwise linear interpolation')
    GL=1;
else
    GL=0;
    GLx=[]; GLy=[] ; xGL=[] ; yGL=[]; xglc=[] ; yglc=[]; 
end

x=coordinates(:,1) ; y=coordinates(:,2);


%% Velocity field

figure
N=1;
CtrlVar.MinSpeedWhenPlottingVelArrows=10;
%CtrlVar.VelColorMap='hot';
QuiverColorGHG(x(1:N:end),y(1:N:end),u(1:N:end),v(1:N:end),CtrlVar)

%%
hold on
%plot(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,'k','LineWidth',2) ; 
%plot(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,'r') ; 
plot(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,'k') ; 
%plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',2);
title(sprintf('t=%-g ',time)) ; xlabel('x (km)') ; ylabel('y (km)')
% 
% % retrograde GL modifications
% slope0a=973.6686e3;  slope0b=1265.712e3;
% hold on ; plot([slope0a/CtrlVar.PlotXYscale,slope0a/CtrlVar.PlotXYscale],[min(y) max(y)]/CtrlVar.PlotXYscale,'m')
% hold on ; plot([slope0b/CtrlVar.PlotXYscale,slope0b/CtrlVar.PlotXYscale],[min(y) max(y)]/CtrlVar.PlotXYscale,'m')
% % 
% axis([900 1600 -120 120])
% figPosition=[50 150 900 500];
% set(gcf,'Position',figPosition)
% text(1500,100,'Ice shelf','FontSize',13,'FontWeight','bold')
% text(1100,100,'Ice sheet','FontSize',13,'FontWeight','bold')
% title(' ')
% %%




%%

if GL
    figure ; PlotFEmesh(coordinates,connectivity,CtrlVar); hold on
    hold on ;plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2);
    hold off
end

%velscale=100000 ; headscale=0.3; sharp=0.3; head=1; col='k'; lw=1;
%hold off;
%figure ;ghg_arrow(x,y,u-mean(u(:)),v,velscale,headscale,sharp,head,col,lw);

%save TestSave
%% 
if nargin >30
figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,a.*rho/1000,'EdgeColor','none')  ; title('Specific Annual Mass Balance  (m w.e.q/a)')
figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,as.*rho/1000,'EdgeColor','none') ; title('Specific Annual Surface Mass Balance')
figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,ab.*rho/1000,'EdgeColor','none') ; title('Specific Annual Basal Mass Balance')
end

%% Surface, Bed, Bedrock
figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,s,'EdgeColor','none') ; %axis equal
hold on ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,b,'EdgeColor','none') ;
hold on ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,B,'EdgeColor','none') ;
lightangle(-45,30) ; lighting phong ;
view(35,20)
hold on ;

%Draftz=Grid1toGrid2(DTxy,s,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
GLz=Grid1toGrid2(DTxy,s,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2)
GLz=Grid1toGrid2(DTxy,B,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2)

title(['surface, bed and bedrock at t=',num2str(time)]) ;xlabel('x') ; ylabel('y') ;colorbar
    
indh0=find(h<=CtrlVar.ThickMin);
if ~isempty(indh0)  % plot locations of <0 thickness
    fprintf(' found %-i thickness values having min thick of %-g or less \n ',numel(indh0),CtrlVar.ThickMin)
    plot3(x(indh0)/CtrlVar.PlotXYscale,y(indh0)/CtrlVar.PlotXYscale,h(indh0),'+r')
end
tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
xlabel('x (km)')
ylabel('y (km)')
zlabel('(m a.s.l.)')
Cbar=colorbar ; title(Cbar,'(m a.s.l.)')
view(-18,20)

%% GF node
if GL
    CtrlVar.PlotXYscale=1000;
    figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,GF.node,'EdgeColor','none') ;  title(' GF.node ') ;lightangle(-45,30);lighting phong ; view(35,20)
    hold on ; plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,1.1+GLx*0,'r')
    %hold on ; plot3(GLxUpper/CtrlVar.PlotXYscale,GLyUpper/CtrlVar.PlotXYscale,1.1+GLxUpper*0,'g','LineWidth',1)
    %hold on ; plot3(GLxLower/CtrlVar.PlotXYscale,GLyLower/CtrlVar.PlotXYscale,1.1+GLxLower*0,'c','LineWidth',1)
end


%%


if GL
    figure
    trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,dGFdt,'EdgeColor','none') ;  title('rate of change in grounding/floating mask ')
    lightangle(-45,30);lighting phong ; view(35,20)
    GLz=Grid1toGrid2(DTxy,dGFdt,GLx,GLy,CtrlVar); hold on ; plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2) ; colorbar
end


%% Surface
figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,s,'EdgeColor','none') ;
lightangle(-45,30);lighting phong ; view(35,20)
xlabel('x') ; ylabel('y') ; colorbar ; %axis equal
hold on ;
%Draftz=Grid1toGrid2(DTxy,s,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
GLz=Grid1toGrid2(DTxy,s,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2)
title(sprintf('s at t=%-g ',time))

if ~isempty(indh0)  
    plot3(x(indh0)/CtrlVar.PlotXYscale,y(indh0)/CtrlVar.PlotXYscale,h(indh0),'+r')
end
tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);

%% Bedrock
figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,B,'EdgeColor','none')  ;
view(35,20) ;lightangle(-45,30) ; lighting phong ;
hold on
%Draftz=Grid1toGrid2(DTxy,B,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
GLz=Grid1toGrid2(DTxy,B,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2)
xlabel('x') ; ylabel('y') ; colorbar ; %axis equal
title(sprintf('B at t=%-g ',time))
tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
 
%% Ice thickness

figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,h,'EdgeColor','none') ;
view(35,20); lightangle(-45,30) ; lighting phong ;
xlabel('x') ; ylabel('y') ; colorbar ; %axis equal
hold on
%Draftz=Grid1toGrid2(DTxy,h,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
GLz=Grid1toGrid2(DTxy,h,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
title(sprintf('h at t=%-g ',time))
tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);


%%

if nargin==10 ; return ; end

%% Basal slipperiness

figure ;
if CtrlVar.CisElementBased
    PlotElementBasedQuantities(connectivity,coordinates/CtrlVar.PlotXYscale,C);
    hold on ; plot(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,'r','LineWidth',2);
else
    %%
    trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,C,'EdgeColor','none') ;
    view(35,20); lightangle(-45,30) ; lighting phong ;
    %%
    hold on
    GLz=Grid1toGrid2(DTxy,C,GLx,GLy,CtrlVar);  plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
end
hold on ;

title('C') ;xlabel('x') ; ylabel('y') ; colorbar ; 
tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
%axis equal
%Draftz=Grid1toGrid2(DTxy,C,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
%% rate factor


figure ;
if CtrlVar.AGlenisElementBased
    PlotElementBasedQuantities(connectivity,coordinates/CtrlVar.PlotXYscale,AGlen);
    hold on ;
    if GL
        plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2);
    end
else
    trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,AGlen,'EdgeColor','none') ;
    view(35,20); lightangle(-45,30) ; lighting phong ;
    hold on
    if GL
        GLz=Grid1toGrid2(DTxy,AGlen,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
    end
end
title('AGlen') ;xlabel('x') ; ylabel('y') ; colorbar ; %axis equal
tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);

%Draftz=Grid1toGrid2(DTxy,AGlen,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)

if ~isempty(dhdt)
    figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,dhdt,'EdgeColor','none') ;
    lightangle(-45,30) ; lighting phong ;
    xlabel('x') ; ylabel('y') ; colorbar ; %axis equal
    title(sprintf('dh/dt at t=%-g ',time))
    if GL
        hold on ; GLz=Grid1toGrid2(DTxy,dhdt,GLx,GLy,CtrlVar);
        plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
    end
    tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
end


figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,u,'EdgeColor','none') ; title(' u ')  ;xlabel('x') ; ylabel('y') ; colorbar ; %axis equal
title(sprintf('u at t=%-g ',time))
view(35,20); lightangle(-45,30) ; lighting phong ;
tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);

figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,v,'EdgeColor','none') ; title(' v ') ; xlabel('x') ; ylabel('y') ; colorbar ; %axis equal
title(sprintf('v at t=%-g ',time))
view(35,20); lightangle(-45,30) ; lighting phong ;
tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);



speed=sqrt(u.*u+v.*v);
figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,speed,'EdgeColor','none') ;
view(35,20); lightangle(-45,30) ; lighting phong ;
title(' speed ') ; xlabel('x') ; ylabel('y') ; colorbar; %axis equal

if GL
    hold on ;
    GLz=Grid1toGrid2(DTxy,speed,GLx,GLy,CtrlVar);
    plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
end
%Draftz=Grid1toGrid2(DTxy,speed,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
title(sprintf('speed at t=%-g ',time))
tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);

%%

figure ; trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,wSurf,'EdgeColor','none') ;
view(35,20); lightangle(-45,30) ; lighting phong ;
if GL
    hold on ; GLz=Grid1toGrid2(DTxy,wSurf,GLx,GLy,CtrlVar);
    plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
end
%Draftz=Grid1toGrid2(DTxy,w,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
title(' wSurf ') ; xlabel('x') ; ylabel('y') ; colorbar; %axis equal
title(sprintf('wSurf at t=%-g ',time))
tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);



% 
% %% Pertubed Velocity field
% N=1;
% 
% u1=u(1) ; v1=v(1);
% MyScale=1;
% figure ;
% quiver(x(1:N:end)/CtrlVar.PlotXYscale,y(1:N:end)/CtrlVar.PlotXYscale,MyScale*(u(1:N:end)-mean(u)),MyScale*(v(1:N:end)-mean(v))) ; axis equal tight
% hold on ; plot(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,'r','LineWidth',2)
% %hold on ; plot(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,'m','LineWidth',1)
% %hold on ; plot(GLxUpper/CtrlVar.PlotXYscale,GLyUpper/CtrlVar.PlotXYscale,'g','LineWidth',1)
% %hold on ; plot(GLxLower/CtrlVar.PlotXYscale,GLyLower/CtrlVar.PlotXYscale,'c','LineWidth',1)
% title(sprintf('Velocity perturbations t=%-g ',time))
% %hold on ; quiver(x(1)/CtrlVar.PlotXYscale,y(1)/CtrlVar.PlotXYscale,MyScale*u(1),MyScale*v(1),'color','r') ; axis equal tight
% u(1)=u1 ; v(1)=v1;
% hold on
% velaxis=axis;
% 





%% plot velocity field onto an image of Antarctica


if CtrlVar.PlotBackgroundImage==1;
    CurDir=pwd; goto_home_directory; cd Antarctica' Images'\
    fprintf(' Loading Lima Composite ')
    load AntarcticLimaComposite.mat
    fprintf(' done \n ')
    cd(CurDir)
    
    figure ; mapshow(AntarticaLimaComposite,R./CtrlVar.PlotXYscale) ;
    axis(velaxis)
    hold on
    quiver(x(1:N:end)/CtrlVar.PlotXYscale,y(1:N:end)/CtrlVar.PlotXYscale,MyScale*u(1:N:end),MyScale*v(1:N:end)) ; axis equal tight
    plot(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,'r','LineWidth',1)
    
end

%% plot fe mesh
% figure
% CPlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar);
% hold on ; plot(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,'r','LineWidth',1)
% hold on ; plot(GLxUpper/CtrlVar.PlotXYscale,GLyUpper/CtrlVar.PlotXYscale,'g','LineWidth',1)
% hold on ; plot(GLxLower/CtrlVar.PlotXYscale,GLyLower/CtrlVar.PlotXYscale,'c','LineWidth',1)
% axis equal tight

%% Plot integration point fields


% in contrast to the nodal coordinates the int coordinates can contain duplicates so I must get rid of them
if CtrlVar.PlotStrains==1;
    
    Xint=xint(:) ; Yint=yint(:);
    [~, Iint, ~] = unique([Xint Yint],'first','rows');
    Iint = sort(Iint); Xint = Xint(Iint); Yint = Yint(Iint);
    DTint = DelaunayTri(Xint,Yint);
    ic=incenters(DTint);
    [cn,on] = inpoly(ic,MeshBoundaryCoordinates);
    DTintTriInside=DTint.Triangulation(cn,:);
    
    etaInt=etaInt(:) ; etaInt = etaInt(Iint); figure ; trisurf(DTintTriInside,Xint/CtrlVar.PlotXYscale,Yint/CtrlVar.PlotXYscale,etaInt,'EdgeColor','none') ; title('eta effective (integration point values)')
    view(35,20); lightangle(-45,30) ; lighting phong ;  hold on ; GLz=Grid1toGrid2(DTint,etaInt,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
    
    exx=exx(:) ; exx = exx(Iint); figure ; trisurf(DTintTriInside,Xint/CtrlVar.PlotXYscale,Yint/CtrlVar.PlotXYscale,exx,'EdgeColor','none') ; title('exx (integration point values)')
    view(35,20); lightangle(-45,30) ; lighting phong ;  hold on ; GLz=Grid1toGrid2(DTint,exx,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
    
    eyy=eyy(:) ; eyy = eyy(Iint); figure ; trisurf(DTintTriInside,Xint/CtrlVar.PlotXYscale,Yint/CtrlVar.PlotXYscale,eyy,'EdgeColor','none') ; title('eyx (integration point values)')
    view(35,20); lightangle(-45,30) ; lighting phong ;  hold on ; GLz=Grid1toGrid2(DTint,eyy,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
    
    exy=exy(:) ; exy = exy(Iint); figure ; trisurf(DTintTriInside,Xint/CtrlVar.PlotXYscale,Yint/CtrlVar.PlotXYscale,exy,'EdgeColor','none') ; title('exy (integration point values)')
    view(35,20); lightangle(-45,30) ; lighting phong ;  hold on ; GLz=Grid1toGrid2(DTint,exy,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
    
    e=e(:) ; e = e(Iint); figure ; trisurf(DTintTriInside,Xint/CtrlVar.PlotXYscale,Yint/CtrlVar.PlotXYscale,e,'EdgeColor','none') ; title('eps effective (integration point values)');
    view(35,20); lightangle(-45,30) ; lighting phong ;  hold on ; GLz=Grid1toGrid2(DTint,e,GLx,GLy,CtrlVar); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
    
    %wint=wint(:) ; wint = wint(Iint); figure ; trisurf(DTintTriInside,Xint/CtrlVar.PlotXYscale,Yint/CtrlVar.PlotXYscale,wint) ; title('w at integration points')
    
    %	wint=wint2(:) ; wint = wint(Iint); figure ; trisurf(DTintTriInside,Xint/CtrlVar.PlotXYscale,Yint/CtrlVar.PlotXYscale,wint) ; title('w2 at integration points')
    
end

%figure ; quiver(x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,u-mean(u),v-mean(v)) ; axis equal tight ; title('Velocity pertubations ')

return





%%
figure ; tricontour(coordinates,TRIxy,b,15) ;  title(' b ') ;xlabel('x') ; ylabel('y');
axis([min(x) max(x) min(y) max(y)]) ; axis equal tight

figure ; tricontour(coordinates,TRIxy,B,15) ;  title(' B ') ;xlabel('x') ; ylabel('y');
axis([min(x) max(x) min(y) max(y)]) ; axis equal tight

figure ; tricontour(coordinates,TRIxy,s,15) ;  title(' s ') ;xlabel('x') ; ylabel('y');
axis([min(x) max(x) min(y) max(y)]) ; axis equal tight

figure ; tricontour(coordinates,TRIxy,h,15) ;  title(' h ') ;xlabel('x') ; ylabel('y');
axis([min(x) max(x) min(y) max(y)]) ; axis equal tight

figure ; tricontour(coordinates,TRIxy,u,15) ;  title(' u ') ;xlabel('x') ; ylabel('y');
axis([min(x) max(x) min(y) max(y)]) ; axis equal tight

figure ; tricontour(coordinates,TRIxy,v,15) ;  title(' v ') ;xlabel('x') ; ylabel('y');
axis([min(x) max(x) min(y) max(y)]) ; axis equal tight

figure ; tricontour(coordinates,TRIxy,w,15) ;  title(' w ') ;xlabel('x') ; ylabel('y');
axis([min(x) max(x) min(y) max(y)]) ; axis equal tight

figure ; tricontour(coordinates,TRIxy,C,15) ;  title(' C ') ;xlabel('x') ; ylabel('y');
axis([min(x) max(x) min(y) max(y)]) ; axis equal tight


qx=linspace(min(x),max(x),100); qy=linspace(min(y),max(y),100);
qz=gridfit(x,y,h,qx,qy);
figure; contour(qx,qy,qz);  title(' h ') ;xlabel('x') ; ylabel('y');


%%



figure ; trisurf(DTintTriInside,xint(:),yint(:),exx(:)) ; title('exx')
figure ; trisurf(DTintTriInside,xint(:),yint(:),eyy(:)) ; title('eyy')
figure ; trisurf(DTintTriInside,xint(:),yint(:),exy(:)) ; title('exy')
e=sqrt(exx.^2+eyy.^2+exx.*eyy+exy.^2);
figure ; trisurf(DTintTriInside,xint(:),yint(:),e(:)) ; title('effective strain rates')

%% not sure why but I could not get tricontour to work with the xint, yint data

TRIint = DelaunayTri(xint(:),yint(:));
[qx,qy] = meshgrid(linspace(min(x),max(x),100),linspace(min(y),max(y),100));

F = TriScatteredInterp(TRIint, e(:),'natural'); qz=F(qx,qy);
F.Method = 'natural';
figure ; [cc,hh]=contourf(qx,qy,qz); title('effective strain rates'); colorbar

F = TriScatteredInterp(TRIint, exx(:),'natural'); qz=F(qx,qy);
figure ; [cc,hh]=contourf(qx,qy,qz); title('exx'); colorbar

F = TriScatteredInterp(TRIint, eyy(:),'natural'); qz=F(qx,qy);
figure ; [cc,hh]=contourf(qx,qy,qz); title('eyy'); colorbar

F = TriScatteredInterp(TRIint, exy(:),'natural'); qz=F(qx,qy);
figure ; [cc,hh]=contourf(qx,qy,qz); title('exy'); colorbar

end
