function FileName=PlotRestartFile(FileName)

persistent LastFileName

if nargin==0
    if isempty(LastFileName)
        [FileName,PathName,FilterIndes]=uigetfile('*.mat','Select a restart file');
    else
        FileName=LastFileName;
    end
end

if isequal(FileName,0) ; return ; end

fprintf(' loading %s ',FileName)
load(FileName)
fprintf(' done \n')

LastFileName=FileName ;




%%


plots='-dhdt-q-h-s-b-';
plots='-ubvb-';
%plots='-sbB-';

CtrlVar=CtrlVarInRestartFile;
coordinates=MUA.coordinates ; connectivity=MUA.connectivity;
GF=GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
tri=TriFE(MUA.connectivity);
speed=sqrt(ub.*ub+vb.*vb);
x=coordinates(:,1); y=coordinates(:,2);
%[GLx,GLy,GLxUpper,GLyUpper,GLxLower,GLyLower] = FindGL(DTxy,GF.node,CtrlVar);
GLgeo=GLgeometry(connectivity,coordinates,GF,CtrlVar);

x=coordinates(:,1); y=coordinates(:,2);


% vel

%%
if ~isempty(strfind(plots,'-ubvb-'))
  
    
    
    figure(10) ; hold off
    N=1;
    
    CtrlVar.VelPlotIntervalSpacing='log10';
    QuiverColorGHG(coordinates(1:N:end,1),coordinates(1:N:end,2),ub(1:N:end),vb(1:N:end),CtrlVar);
    
    
    if ~isempty(strfind(plots,'-ShowBoundaryNodes-'))
        [Boundary.Nodes,Boundary.EdgeCornerNodes,Boundary.FreeElements,Boundary.Edges]=FindBoundaryNodes(connectivity,coordinates);
        plot(coordinates(Boundary.Nodes,1)/CtrlVar.PlotXYscale,coordinates(Boundary.Nodes,2)/CtrlVar.PlotXYscale,'.r')
    end
    
    title(sprintf('(u,v) t=%-g ',time)) ; xlabel('x (km)') ; ylabel('y (km)') ; axis equal tight
    %%
    hold off
    
end

%%
if ~isempty(strfind(plots,'-q-'))
    figure(20) ; hold off
    N=1;
    %CtrlVar.MinSpeedWhenPlottingVelArrows=1; CtrlVar.MaxPlottedSpeed=3000; CtrlVar.VelPlotIntervalSpacing='log10';
    %CtrlVar.VelColorMap='hot';
    
    qx=u.*h;  qy=v.*h;
    PlotFEmesh(coordinates,connectivity,CtrlVar) ; hold on
    CtrlVar.VelColorBarTitle='(m^2 a^{-1})' ;
    QuiverColorGHG(coordinates(1:N:end,1),coordinates(1:N:end,2),qx(1:N:end),qy(1:N:end),CtrlVar);
    hold on
    
    title(sprintf('q at t=%-g ',time)) ; xlabel('x (km)') ; ylabel('y (km)')
    
    hold off
    
end

if ~isempty(strfind(plots,'-qgrad-')) ||  ~isempty(strfind(plots,' dhdtest '))
    
    
    
    qx=u.*h;  qy=v.*h;
    [dqxdx,dqxdy,xint,yint]=calcFEderivatives(qx,coordinates,connectivity,nip,CtrlVar);
    [dqydx,dqydy,xint,yint]=calcFEderivatives(qx,coordinates,connectivity,nip,CtrlVar);
    
    [DTxy,tri,DTint,DTintTriInside,Xint,Yint,xint,yint,Iint]=TriangulationNodesIntegrationPoints(coordinates,connectivity,Boundary.EdgeCornerNodes,nip);
    
    temp=dqxdx(:) ; temp=temp(Iint); dqxdxNode=Grid1toGrid2(DTint,temp,x,y);
    temp=dqydy(:) ; temp=temp(Iint); dqydyNode=Grid1toGrid2(DTint,temp,x,y);
    
    
    if ~isempty(strfind(plots,' qgrad '))
        figure(30) ; hold off
        
        PlotFEmesh(coordinates,connectivity,CtrlVar)
        hold on
        CtrlVar.VelColorBarTitle='(m a^{-1})' ;
        
        QuiverColorGHG(xint(:),yint(:),dqxdx(:),dqydy(:),CtrlVar);
        %QuiverColorGHG(x,y,dqxdxNode,dqydyNode,CtrlVar);
        hold on
        plot(x(Boundary.EdgeCornerNodes)/CtrlVar.PlotXYscale,y(Boundary.EdgeCornerNodes)/CtrlVar.PlotXYscale,'k')
        
        title(sprintf('grad q at t=%-g ',time)) ; xlabel('x (km)') ; ylabel('y (km)')
        I=find(h<=CtrlVar.ThickMin);
        plot(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,'.r')
        hold off
        
    end
    
    if   ~isempty(strfind(plots,'-dhdtest-'))
        
        
        
        
        
        dhdtEst=a-dqxdxNode-dqydyNode;
        
        I=find(h<=CtrlVar.ThickMin);  dhdtEst(I)=dhdt(I);
        
        figure(40) ; hold off
        trisurf(tri,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,dhdtEst,'EdgeColor','none') ;
        view(0,90); lightangle(-45,30) ; lighting phong ; axis equal tight
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        colorbar ; title(colorbar,'(m)')
        title(sprintf('estimated dhdt at t=%#5.1f ',time))
        hold off
        
        figure(50) ; hold off
        trisurf(tri,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,dhdtEst-dhdt,'EdgeColor','none') ;
        view(0,90); lightangle(-45,30) ; lighting phong ; axis equal tight
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        colorbar ; title(colorbar,'(m)')
        title(sprintf('error in dhdt at t=%#5.1f ',time))
        hold off
        
        fprintf(' std(dhdtEst-dhdt)=%-g \t min(dhdtEst-dhdt)=%-g \t max(dhdtEst-dhdt)=%-g \n',std(dhdtEst-dhdt),min(dhdtEst-dhdt),max(dhdtEst-dhdt))
    end
end


% ice thickness

if ~isempty(strfind(plots,'-speed-'))
    %%
    speed=sqrt(u.*u+v.*v);
    figure(65) ; hold off ; PlotNodalBasedQuantities(tri,coordinates,speed,CtrlVar);
    xlabel('(km)') ; ylabel('(km)') ;
    title(colorbar,'(m)') ; %caxis([0 100])
    hold on
    
    title(sprintf('speed at t=%#5.1f ',time))
    %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
    axis equal tight
    
end

if ~isempty(strfind(plots,'-h-'))
    %%
    
    figure(60) ; hold off ; 
    PlotMeshScalarVariable(CtrlVar,MUA,h);
    xlabel('xps (km)') ; ylabel('yps (km)') ;
    title(colorbar,'(m)') ; %caxis([0 100])
    hold on
    
    title(sprintf('h at t=%#5.1f ',time))
    %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
    axis equal tight
    if ~isempty(strfind(plots,'-ShowMinThick-'))
        I=find(h<=CtrlVar.ThickMin);
        plot3(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,h(I)+10,'.r')
    end
    hold off
    %%
end

%%
if ~isempty(strfind(plots,'-sbB-'))
    %%
    
    hold off
    trisurf(tri,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,B,'EdgeColor','none') ;
    hold on
    trisurf(tri,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,s,'EdgeColor','none') ;
    trisurf(tri,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,b,'EdgeColor','none') ;
    
    view(30,90); lightangle(-45,30) ; lighting phong ;
    xlabel('x (km)') ; ylabel('y(km)') ;
    colorbar ; title(colorbar,'(m)')

end
%%


if ~isempty(strfind(plots,'-as-'))
    %%
    
    figure(160) ; hold off ; 
    PlotMeshScalarVariable(CtrlVar,MUA,as);
    xlabel('(km)') ; ylabel('(km)') ;
    title(colorbar,'(m w.eq.)') ; %caxis([0 100])
    hold on
    
    title(sprintf('as at t=%#5.1f ',time))
    %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
    axis equal tight
    
    %%
end

if ~isempty(strfind(plots,'-s-'))
    %  surface
    figure(70) ; hold off
    PlotMeshScalarVariable(CtrlVar,MUA,s);
    
    xlabel('(km)') ; ylabel('(km)') ;
    colorbar ; title(colorbar,'(m)')
    hold on
    
    title(sprintf('s at t=%#5.1f ',time))
    %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
    axis equal tight
    hold off
    
end

if ~isempty(strfind(plots,'-a-'))
    %  surface
    figure(71) ; hold off
    
    PlotMeshScalarVariable(CtrlVar,MUA,a);
    xlabel('(km)') ; ylabel('(km)') ;
    colorbar ; title(colorbar,'(m/a)')
    hold on
    
    title(sprintf('a at t=%#5.1f ',time))
    %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
    axis equal tight
    hold off
    
end

if ~isempty(strfind(plots,'-B-'))
    %  surface
    figure(80) ; hold off
    trisurf(tri,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,B,'EdgeColor','none') ;
    view(0,90); lightangle(-45,30) ; lighting phong ;
    xlabel('xps (km)') ; ylabel('yps (km)') ;
    colorbar ; title(colorbar,'(m)')
    hold on
    
    title(sprintf('B at t=%#5.1f ',time))
    %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
    axis equal tight
    hold off
    
end

% dh/dt
if ~isempty(strfind(plots,'-dhdt-'))
    figure(90) ; hold off
    trisurf(tri,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,dhdt,'EdgeColor','none') ;
    view(0,90); lightangle(-45,30) ; lighting phong ;
    xlabel('xps (km)') ; ylabel('yps (km)') ;
    colorbar ; title(colorbar,'(m)') ;
    hold on
    
    title(sprintf('dh/dt at t=%#5.1f ',time))
    %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
    axis equal tight
    
    if ~isempty(strfind(plots,'-ShowMinThick-'))
        I=find(h<=CtrlVar.ThickMin);
        plot3(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,h(I)+10,'.r')
    end
    hold off
end

if ~isempty(strfind(plots,'-boundary-'))
    
    PlotBoundary(Boundary,connectivity,coordinates,CtrlVar);
    
end




end


