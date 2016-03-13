function PlotResultsFile(FileName,plots)
    
    %% PlotResultsFile(FileName,plots)
    %  Plots data in `FileName'
    %  plots  : a sting indicating which fields to plot
    %           for example: string='-uv-h-' plots the velocity vectors (u,v) and ice thickness h
    %
    %
    
    
    if nargin==0
        [FileName,PathName,FilterIndes]=uigetfile('*.mat');
        
    end
    
    if nargin<2
        plots='-uv-' ;
    end
    
    if isequal(FileName,0) ; return ; end
    
    fprintf(' Plotting results in file %-s \n ',FileName)
    load(FileName,'CtrlVar','TRIxy','coordinates','connectivity','s','h','b','B','u','v','wSurf','dhdt','dsdt','time','GF','rho','rhow','as','ab')
    
    
     [Boundary.Nodes,Boundary.EdgeCornerNodes,Boundary.FreeElements,Boundary.Edges,Boundary.Edge]=FindBoundaryNodes(connectivity,coordinates);
     GLgeo=GLgeometry(connectivity,coordinates,GF,CtrlVar);
     x=coordinates(:,1);  y=coordinates(:,2); 
     speed=sqrt(u.*u+v.*v);
     tri=TriFE(connectivity);
    
     %% -ab-  basal melt

     if ~isempty(strfind(plots,'-ab-'))
         figure
         
         [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,coordinates,ab,CtrlVar);
         xlabel('xps (km)') ; ylabel('yps (km)') ;
         colorbar ; title(colorbar,'(m/a)')
         hold on
         
         plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
         title(sprintf('a_b at t=%#5.1f ',time))
         axis equal tight
         hold off
     end
     
     
     
     
     %% -sbB-
    
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
        
        title(sprintf('t=%#5.1f ',time))
        axis equal
        tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) 1000*tt(3)]);
        axis tight
        hold off
    end
    
    %% -uv-
    
    if ~isempty(strfind(plots,'-uv-'))
        figure(10)
        N=1;
        CtrlVar.MinSpeedWhenPlottingVelArrows=1; CtrlVar.MaxPlottedSpeed=max(speed); %CtrlVar.VelPlotIntervalSpacing='log10';
        %CtrlVar.VelColorMap='hot';
        QuiverColorGHG(x(1:N:end),y(1:N:end),u(1:N:end),v(1:N:end),CtrlVar);
        hold on
        %plot(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,'k','LineWidth',2) ;
        %plot(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,'r') ;
        %plot(GLx/CtrlVar.PlotXYscale,GLy/CtrlVar.PlotXYscale,'k') ;
        plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2);
        plot(x(Boundary.EdgeCornerNodes)/CtrlVar.PlotXYscale,y(Boundary.EdgeCornerNodes)/CtrlVar.PlotXYscale,'k')
        
        title(sprintf('t=%-g ',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
        axis equal tight
    end
    
    % -speed-
    if ~isempty(strfind(plots,'-speed-'))
        figure(20)
        [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,coordinates,log10(speed),CtrlVar);
                xlabel('xps (km)') ; ylabel('yps (km)') ;
        title(sprintf('log10(speed) at t=%#5.1f ',time))
        colorbar ; ; title(colorbar,'log_{10}(m/a)')
        hold on
        
        plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
        
        axis equal tight
        hold off
    end
    
    
    
    
    %% -C-
    if ~isempty(strfind(plots,'-C-'))
        figure(25)
        hold off
        [hPatch]=PlotElementBasedQuantities(coordinates/CtrlVar.PlotXYscale,connectivity,log10(C));
        
        
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        title(sprintf('log10(Slipperiness) at t=%#5.1f ',time))
        colorbar  ; title(colorbar,'log_{10}(C) (m/(kPa^3 a))')
        hold on
        
        plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
        
        hold off
    end
    
    
    %% -h-  ice thickness
    if ~isempty(strfind(plots,'-h-'))
        figure(30)
        hold off
        [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,coordinates,h,CtrlVar);
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        colorbar ; title(colorbar,'(m)')
        hold on
        
        plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
        title(sprintf('h at t=%#5.1f ',time))
        
        hold off
    end
    
    %% -mesh-
    if ~isempty(strfind(plots,'-mesh-'))
        
        figure(50) ; hold off
        PlotFEmesh(coordinates,connectivity,CtrlVar) ;
        hold on
        plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2);
        xlabel('xps (km)') ; ylabel('yps (km)') ;
    end
    
    
    %% -AGlen-
    if ~isempty(strfind(plots,'-AGlen-'))
        figure(55)
        hold off
        [hPatch]=PlotElementBasedQuantities(coordinates/CtrlVar.PlotXYscale,connectivity,log10(AGlen));
        
        
        xlabel('xps (km)') ; ylabel('yps (km)') ;
        title(sprintf('log10(A) at t=%#5.1f ',time))
        colorbar  ; title(colorbar,'log_{10}(A) (kPa^{-3} a^{-1})')
        hold on
        
        plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
        
        hold off
    end
    

    
    
    
end