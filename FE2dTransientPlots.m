function FE2dTransientPlots(CtrlVar,DTxy,TRIxy,MeshBoundaryCoordinates,GF,dGFdt,coordinates,connectivity,...
        b,B,S,s,h,u,v,wSurf,dhdt,dsdt,dbdt,C,AGlen,m,n,xint,yint,wSurfInt,etaInt,exx,eyy,exy,e,time,...
        rho,rhow,a,as,ab,Boundary,nip)
    

    
    %%
    %plots=' q uv dhdt h qgrad ';
    %plots=' save only ';
    plots='-h(x)-dh/dt(x)-';
    plots='-sbB-uv-speed-GF-';
    
    fprintf('\n --------- FE2dTransientPlots at t=%-g \n \n ',time)
    
    if ~isempty(strfind(plots,'-save-'))
        
        if strcmp(CtrlVar.FE2dTransientPlotsInfostring,'Last transient plot')==0
            
            FileName=['ResultsFiles/',sprintf('%07i',round(100*time)),'-TransPlots-',CtrlVar.Experiment];
            fprintf(' Saving data in %s \n',FileName)
            save(FileName,'CtrlVar','TRIxy','coordinates','connectivity','s','h','b','B','u','v','wSurf','dhdt','dsdt','time','GF','rho','rhow','as','ab')
            
        end
    end
    
    if ~isempty(strfind(plots,' save only '))
        return
    end
    
    speed=sqrt(u.*u+v.*v);
    x=coordinates(:,1); y=coordinates(:,2);
    %[GLx,GLy,GLxUpper,GLyUpper,GLxLower,GLyLower] = FindGL(DTxy,GF.node,CtrlVar);
    
    
    GLgeo=GLgeometry(connectivity,coordinates,GF,CtrlVar);
    
    %% h(x) 
    if ~isempty(strfind(plots,'-h(x)-'))
        %  lower surface
        figure(140)
        hold off
        x=coordinates(:,1) ;
        [xsorted,I]=sort(x);
        plot(x(I)/CtrlVar.PlotXYscale,h(I),'.-')
        xlabel('xps (km)') ; ylabel('h (m)') ;
        hold on
        title(sprintf('h(x) at t=%#5.1g ',time))
        hold off
        
    end
    
      %% dh/dt(x) , good for
    if ~isempty(strfind(plots,'-dh/dt(x)-'))
        %  lower surface
        figure(150)
        hold off
        x=coordinates(:,1) ;
        [xsorted,I]=sort(x);
        plot(x(I)/CtrlVar.PlotXYscale,dhdt(I),'.-')
        xlabel('xps (km)') ; ylabel('dh/dt (m/yr)') ;
        hold on
        title(sprintf('dhdt(x) at t=%#5.1g ',time))
        hold off
        
    end
    
    if ~isempty(strfind(plots,'-GF-'))
        %  lower surface
        figure(40)
        hold off
        PlotNodalBasedQuantities(connectivity,coordinates,GF.node,CtrlVar);
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        colorbar ; title(colorbar,'(m)')
        hold on
        
        plot3(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,0*GLgeo(:,[3 4])'+10000,'r','LineWidth',2)
        %Draftz=Grid1toGrid2(DTxy,h,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
        %GLz=Grid1toGrid2(DTxy,h,GLx,GLy); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
        title(sprintf('GF.node at t=%#5.1g ',time))
        %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
        axis equal tight
        hold off
        
    end
    
    
    if ~isempty(strfind(plots,'-s-'))
        %  lower surface
        figure(40)
        hold off
        PlotNodalBasedQuantities(connectivity,coordinates,s,CtrlVar);
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        colorbar ; title(colorbar,'(m)')
        hold on
        
        plot3(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,0*GLgeo(:,[3 4])'+10000,'r','LineWidth',2)
        %Draftz=Grid1toGrid2(DTxy,h,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
        %GLz=Grid1toGrid2(DTxy,h,GLx,GLy); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
        title(sprintf('s at t=%#5.1g ',time))
        %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
        axis equal tight
        hold off
        
    end
    %% sbB
    
    if ~isempty(strfind(plots,'-sbB-'))
        figure(5)
        hold off
        trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,s,'EdgeColor','none') ; hold on
        trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,b,'EdgeColor','none') ;
        trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,B,'EdgeColor','none') ;
        view(50,20); lightangle(-45,30) ; lighting phong ;
        xlabel('y') ; ylabel('x') ;
        colorbar ; title(colorbar,'(m)')
        hold on
        
        title(sprintf('sbB at t=%#5.1g ',time))
        axis equal ; tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) 1000*tt(3)/CtrlVar.PlotXYscale]); axis tight
        hold off
    end
    
    %%
%     plots='-uv-';
%     speed=sqrt(u.*u+v.*v);
%     x=coordinates(:,1); y=coordinates(:,2);
%     GF=GL2d(B,S,h,rhow,rho,connectivity,CtrlVar);
%     GLgeo=GLgeometry(connectivity,coordinates,GF,CtrlVar);
%     [Boundary.Nodes,Boundary.EdgeCornerNodes]=FindBoundaryNodes(connectivity,coordinates);
    if ~isempty(strfind(plots,'-uv-'))
        figure(11)
        N=1;
        CtrlVar.MinSpeedWhenPlottingVelArrows=0; CtrlVar.MaxPlottedSpeed=max(speed); 
        CtrlVar.VelPlotIntervalSpacing='log10';
        %CtrlVar.VelColorMap='hot';
        QuiverColorGHG(x(1:N:end),y(1:N:end),u(1:N:end),v(1:N:end),CtrlVar);
        hold on
        %plot(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,'k','LineWidth',2) ;
        %plot(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,'r') ;
        %plot(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,'k') ;
        plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2);
        plot(x(Boundary.EdgeCornerNodes)/CtrlVar.PlotXYscale,y(Boundary.EdgeCornerNodes)/CtrlVar.PlotXYscale,'k')
        
        title(sprintf('(u,v) t=%-g ',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
        axis equal tight
    end
    
    % speed
    if ~isempty(strfind(plots,'-speed-'))
        figure(20)
        trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,log10(speed),'EdgeColor','none') ;
        view(0,90); lightangle(-45,30) ; lighting phong ;
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        title(sprintf('log10(speed) at t=%#5.1f ',time))
        colorbar ; caxis([0 4]) ; title(colorbar,'log_{10}(m/a)')
        hold on
        
        plot3(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,0*GLgeo(:,[3 4])'+10000,'r','LineWidth',2)
        %Draftz=Grid1toGrid2(DTxy,h,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
        %GLz=Grid1toGrid2(DTxy,h,GLx,GLy); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
        
        %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
        axis equal tight
        hold off
    end
    
    
     % basal melt rate
    if ~isempty(strfind(plots,'-ab-'))
        figure(25)
        trisurf(TRIxy,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,ab,'EdgeColor','none') ;
        view(0,90); lightangle(-45,30) ; lighting phong ;
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        title(sprintf('Basal melting at t=%#5.1f ',time))
        colorbar  ; title(colorbar,'a_b (m/a)')
        hold on
        
        plot3(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,0*GLgeo(:,[3 4])'+10000,'r','LineWidth',2)
        %Draftz=Grid1toGrid2(DTxy,h,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
        %GLz=Grid1toGrid2(DTxy,h,GLx,GLy); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
        
        %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
        axis equal tight
        hold off
    end
    
      % basal melt rate
    if ~isempty(strfind(plots,'-as-'))
        figure(25)
        PlotNodalBasedQuantities(connectivity,coordinates,as,CtrlVar);
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        title(sprintf('Surface mass balance at t=%#5.1f ',time))
        colorbar  ; title(colorbar,'a_s (m/a)')
        hold on
        
        plot3(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,0*GLgeo(:,[3 4])'+10000,'r','LineWidth',2)
        %Draftz=Grid1toGrid2(DTxy,h,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
        %GLz=Grid1toGrid2(DTxy,h,GLx,GLy); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
        
        %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
        axis equal tight
        hold off
    end
    
    
    
     % C
    if ~isempty(strfind(plots,'-C-'))
        figure(25)
        hold off
        [hPatch]=PlotElementBasedQuantities(coordinates/CtrlVar.PlotXYscale,connectivity,log10(C));
        
        
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        title(sprintf('log10(Slipperiness) at t=%#5.1f ',time))
        colorbar  ; title(colorbar,'log_{10}(C) (m/(kPa^3 a))')
        hold on
        
        plot3(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,0*GLgeo(:,[3 4])'+10000,'r','LineWidth',2)
        %Draftz=Grid1toGrid2(DTxy,h,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
        %GLz=Grid1toGrid2(DTxy,h,GLx,GLy); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
        
        %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
        axis equal tight
        hold off
    end
    
    
    % ice thickness
    if ~isempty(strfind(plots,'-h-'))
        figure(30)
        hold off
        PlotNodalBasedQuantities(connectivity,coordinates,h,CtrlVar);
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        colorbar ; title(colorbar,'(m)')
        hold on
        
        plot3(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,0*GLgeo(:,[3 4])'+10000,'r','LineWidth',2)
        %Draftz=Grid1toGrid2(DTxy,h,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
        %GLz=Grid1toGrid2(DTxy,h,GLx,GLy); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
        title(sprintf('h at t=%#5.1f ',time))
        %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
        axis equal tight
        hold off
    end
    
    if ~isempty(strfind(plots,'-b-'))
        %  lower surface
        figure(40)
        hold off
        PlotNodalBasedQuantities(connectivity,coordinates,b,CtrlVar);
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        colorbar ; title(colorbar,'(m)')
        hold on
        
        plot3(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,0*GLgeo(:,[3 4])'+10000,'r','LineWidth',2)
        %Draftz=Grid1toGrid2(DTxy,h,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
        %GLz=Grid1toGrid2(DTxy,h,GLx,GLy); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
        title(sprintf('b at t=%#5.1f ',time))
        %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
        axis equal tight
        hold off
        
    end
    
    if ~isempty(strfind(plots,'-mesh-'))
        CtrlVar.PlotMesh=1;
        figure(50) ; hold off
        PlotFEmesh(coordinates,connectivity,CtrlVar) ; 
        hold on
        plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2);
        xlabel('xps (km)') ; ylabel('yps (km)') ;
    end
    
    
      % AGlen
    if ~isempty(strfind(plots,'-AGlen-'))
        figure(55)
        hold off
        [hPatch]=PlotElementBasedQuantities(coordinates/CtrlVar.PlotXYscale,connectivity,log10(AGlen));
        
       
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        title(sprintf('log10(A) at t=%#5.1f ',time))
        colorbar  ; title(colorbar,'log_{10}(A) (kPa^{-3} a^{-1})')
        hold on
        
        plot3(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,0*GLgeo(:,[3 4])'+10000,'r','LineWidth',2)
        %Draftz=Grid1toGrid2(DTxy,h,Draftx,Drafty); plot3(Draftx/CtrlVar.PlotXYscale,Drafty/CtrlVar.PlotXYscale,Draftz,'m','LineWidth',2)
        %GLz=Grid1toGrid2(DTxy,h,GLx,GLy); plot3(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,GLz,'r','LineWidth',2);
        
        %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
        axis equal tight
        hold off
    end
    
    
    
    %%
end  
 