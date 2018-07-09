% a simple plotting routine to plot a sequence of results filesm
%% variables:
CreateVideo=1;
pos=[40 200 2300 1800]; % hig res
pos=[200 50 1200 900];

cd C:\cygwin64\home\Hilmar\GHG\Ua\Runs\PIG\PIG_Thwaites
CurDir=pwd;
RunTable=readtable('RunTable.csv');

%Run=RunTable.RunId;

clear Run
Run{1}=RunTable.RunId{2};
%Run{2}=RunTable.RunId{2};
%Run{3}=RunTable.RunId{3};





cd(DriveName2DriveLetter('WorkSSD1')) ; cd GHG\Ua2D-ResultsFiles\PIG-Thwaites


% 
% iRun=1;  Run{iRun}='J2-Tw-100-tcDe700-Sw3460-sin-deltaTw0-PeriodTw100-ahFeedback0Edge-Wise-supg';





plots='-ubvb-';
%plots='-h-';
%plots='-s-';
%plots='-dhdt-';
plots='-ab-';%plots='-mesh-';%plots='-ab-h-';
%plots='-MeltNodes-';
%plots='-log10(BasalSpeed)-';
%plots='-sbB-';
plots='-Bab-';
dt=5;
PlotMinThickLocations=0;   % nodes at min thickness shown as red dots

%%
PlotArea=NaN;
PlotRegion='pigiceshelf';
PlotRegion='pig-twg';
%PlotRegion=[];




for J=1:numel(Run)
    
    usrstr=[]; TRI=[]; DT=[] ; Par=[]; iFile=1; iFrame=0; k=0; 
    
    %   FileName=F{J};
    
    fprintf(' %i \t %s: \n',J,Run{J})
    
    if CreateVideo
        
        %hFig = figure('Visible','off');
        %set(hFig, 'PaperPositionMode','auto')
        %vidObj = VideoWriter([Run{J},plots,'.avi'],'Archival');
        
        
        vidObj = VideoWriter([Run{J},plots,PlotRegion,'.avi']);
        vidObj.FrameRate=2;
        %vidObj.FrameRate=1;   % frames per sec
        open(vidObj);
    end
    
    list=dir(['*',Run{J},'*.mat']);
    nFiles=length(list);
    
    
    %FigHandleh=figure('Position',[100 50 1100 800],'units','pixel') ;
    %FigHandleh=figure('Position',[100 50 1000 700],nits','pixel') ;
    
    
    
   
    %nFrames=100 ; Frame(nFrames) = struct('cdata',[],'colormap',[]);
    
    
    
    
    while iFile<=nFiles
        
        %for I=1:100
        
        %if strcmp(list(I).name(6:7),'00')
        t=str2double(list(iFile).name(1:7))/100;
        if mod(t,dt)==0 && t<=500
            
            
            try
                load(list(iFile).name)
                fprintf(' %s \n ',list(iFile).name)
            catch
                fprintf('could not load %s \n ',list(iFile).name)
            end
            
            if ~isempty(strfind(CtrlVar.Experiment,'Jenkins')) || ~isempty(strfind(CtrlVar.Experiment,'J2'))
                [Tw,Sw,tcDe]=TwSwtcDe(CtrlVar,time);
                fprintf('time=%f \t Tw=%f \t De=%f \n',time,Tw,tcDe)
            else
                Tw=[];
            end
            
            
            GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
            TRI=[]; DT=[]; xGL=[];yGL=[];
            x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2);
            ih=h<=CtrlVar.ThickMin;
            
                        %             if 

                        

            %                 set(gca,'nextplot','replacechildren');
            %                 set(gcf,'Renderer','zbuffer');
            %             end
            %
            if ~isempty(strfind(plots,'-mesh-'))
                hold off
                PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
                title(sprintf('t=%-g (yr)  #Ele=%-i, #Nodes=%-i, #nod=%-i',time,MUA.Nele,MUA.Nnodes,MUA.nod))
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
                
            end
            
            if ~isempty(strfind(plots,'-ubvb-'))
                % plotting horizontal velocities
                %%
                if ~exist('fubvb','var') || ~ishandle(fubvb)
                    fubvb=figure;
                else
                    figure(fubvb)
                end
                
                if CreateVideo
                    fubvb.Position=pos;
                end
                
                
                hold off
                N=1;
                %speed=sqrt(ub.*ub+vb.*vb);
                %CtrlVar.MinSpeedWhenPlottingVelArrows=0; CtrlVar.MaxPlottedSpeed=max(speed); %
                CtrlVar.VelPlotIntervalSpacing='log10';
                %CtrlVar.VelColorMap='hot';
                %CtrlVar.RelativeVelArrowSize=10;
                PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,'b')
                hold on
                QuiverColorGHG(x(1:N:end),y(1:N:end),ub(1:N:end),vb(1:N:end),CtrlVar);
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                title(sprintf('(ub,vb) t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                %%
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
            end
            
            
            
            
            if ~isempty(strfind(plots,'-log10(BasalSpeed)-'))
                %%
                %us=ub+ud;  vs=vb+vd;
                SurfSpeed=sqrt(ub.*ub+vb.*vb);
                hold off
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(SurfSpeed),CtrlVar);
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                
                title(sprintf('log_{10}(Basal speed) t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                title(colorbar,'log_{10}(m/yr)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                %%
            end
            
            if ~isempty(strfind(plots,'-ab-'))
                if ~exist('fab','var') || ~ishandle(fab)
                    fab=figure;
                    fab.NextPlot='replacechildren';
                else
                    close(fab)
                    fab=figure;
                    
                end
                
                %set(gcf,'nextplot','replacechildren')
                %set(gca,'nextplot','replace')
                %fab.NextPlot='replacechildren';
                %fab.NextPlot
                
                fig=gcf;
                ax=gca;
                
                if CreateVideo
                    fab.Position=pos;
                    %fab.NextPlot='replacechildren';  % slows things down
                end
                
                %ax.NextPlot='replace';
                
                ab(GF.node>0.5)=NaN;
                
                hold off
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,ab,CtrlVar);
                
                hold on
                
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                hold on
                PlotMuaBoundary(CtrlVar,MUA,'b')
                hold on
                title(sprintf('Basal melt at t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                title(colorbar,'(m/yr)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
                caxis([min(ab) max(ab)])
                
            end
            
            if ~isempty(strfind(plots,'-Bab-'))
                if ~exist('fab','var') || ~ishandle(fab)
                    fab=figure;
                    fab.NextPlot='replacechildren';
                else
                    close(fab)
                    fab=figure;
                    
                end
                
                
                %set(gcf,'nextplot','replacechildren')
                %set(gca,'nextplot','replace')
                %fab.NextPlot='replacechildren';
                %fab.NextPlot
                
                fig=gcf;
                ax=gca;
                
                fab.Position=pos;
                
                
                %ax.NextPlot='replace';
                
                ab(GF.node>0.5)=NaN;
                
                hold off
                
                subplot(2,2,1,'Position',[0.05 0.55 0.43 0.43])
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,-ab,CtrlVar);
                
                hold on
                
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                hold on
                PlotMuaBoundary(CtrlVar,MUA,'b')
                hold on
                title(sprintf('Basal melt at t=%-g (yr)',time)) ;
                xlabel('xps (km)') ; ylabel('yps (km)')
                title(colorbar,'(m/yr)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
                caxis([0 100])
                PlotLatLonGrid(1000,2,10);
                ax1=gca;
                
                
                ax2=axes('Position',[0.34 0.8 0.04 0.1]);
                ax2.XLim=[0 1];
                
                TwMax=CtrlVar.Tw+CtrlVar.DeltaTw;
                TwMin=CtrlVar.Tw-CtrlVar.DeltaTw;
                
                ax2.YLim=[TwMin-0.5 TwMax+0.5];
                ax2.XAxis.Visible='off';
                ax2.XAxis.LineWidth=2;
                line([0 1],[Tw Tw],'Parent',ax2,'color','r','LineWidth',2)
                ax2.YAxis.Label.String='T_{CDW} (C^\circ)';
                ax2.Color=[0.95 0.95 0.95];
                
                
                
                ax3=axes('Position',[0.34 0.68 0.04 0.1]);
                ax3.XLim=[0 1];
                
               
                temp1=CtrlVar.tcDe+CtrlVar.DeltaDe;
                temp2=CtrlVar.tcDe-CtrlVar.DeltaDe;
                DeMin=min([temp1 ; temp2]);
                DeMax=max([temp1 ; temp2]);
                
                
                ax3.YLim=[DeMin-50 DeMax+50];
                
                
                
                ax3.XAxis.Visible='off';
                ax3.XAxis.LineWidth=2;
                line([0 1],[tcDe tcDe],'Parent',ax3,'color','r','LineWidth',2)
                ax3.YAxis.Label.String='De (m)';
                ax3.Color=[0.95 0.95 0.95];
                
                
                
                subplot(2,2,2,'Position',[0.55 0.55 0.43 0.43])
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,B,CtrlVar);
                hold on
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r','LineWidth',2);
                PlotMuaBoundary(CtrlVar,MUA,'b');
                caxis([-1200 100])
                title(colorbar,'(m)')
                title(sprintf('Bedrock')) ;
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
                xlabel('xps (km)') ; ylabel('yps (km)')
                PlotLatLonGrid(1000);
                hold off
                
                subplot(2,2,3,'Position',[0.05 0.05 0.43 0.43])
                N=10;
                PlotMuaBoundary(CtrlVar,MUA,'b') ; hold on
                
                if isempty(Par)
                    Par.MaxPlottedSpeed=4000;
                    Par.MinPlottedSpeed=100;
                    Par.VelPlotIntervalSpacing='log10';
                    Par.VelColorBarTitle='(m a^{-1})' ;
                    Par.PlotXYscale=1000;
                end
                
                
                [cbar,Qhandle,Par]=QuiverColorGHG(x(1:N:end),y(1:N:end),ub(1:N:end),vb(1:N:end),Par);
                
                hold on
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r','LineWidth',2);
                title(sprintf('Ice velocity at t=%-g (yr)',time)) ;
                SetRegionalPlotAxis(PlotRegion);
                xlabel('xps (km)') ; ylabel('yps (km)')
                
                PlotLatLonGrid(1000);
                
                
                cbar=colorbar;
                cbar.TickLabels=Par.ticklable;
                cbar.Ticks=Par.tickpos*(cbar.Limits(2)-cbar.Limits(1))+cbar.Limits(1);
                
                
                

                
                subplot(2,2,4,'Position',[0.55 0.05 0.43 0.43])
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,h,CtrlVar);
                
                hold on
                
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                hold on
                PlotMuaBoundary(CtrlVar,MUA,'b')
                hold on
                title(sprintf('ice thickness t=%-g (yr)',time)) ;
                
                title(colorbar,'(m)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
                caxis([100 800])
                xlabel('xps (km)') ; ylabel('yps (km)')
                PlotLatLonGrid(1000);
                
                %set(gcf,'NextPlot','add');
                %axes;
                %text(0.5,0.98,sprintf('time=%-g',time),'Units','normalized','color','r')
                %set(gca,'Visible','off');
                
                
                
                
            end
            
            if ~isempty(strfind(plots,'-b-'))
                if ~exist('fab','var') || ~ishandle(fab)
                    fab=figure;
                else
                    figure(fab)
                end
                if CreateVideo
                    fab.Position=pos;
                end
                hold off
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,b,CtrlVar);
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                
                title(sprintf('b at t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                title(colorbar,'(m)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
                caxis([-400 0])
            end
            
            if ~isempty(strfind(plots,'-s-'))
                if ~exist('fas','var') || ~ishandle(fas)
                    fas=figure;
                else
                    figure(fas)
                end
                if CreateVideo
                    fab.Position=pos;
                end
                hold off
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,s,CtrlVar);
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                
                title(sprintf('s at t=%-g (yr)',time)) ; 
                xlabel('xps (km)') ; ylabel('yps (km)')
                title(colorbar,'(m)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
                %caxis([0 500])
            end
            
            if ~isempty(strfind(plots,'-dhdt-'))
                if ~exist('fas','var') || ~ishandle(fas)
                    fas=figure;
                else
                    figure(fas)
                end
                if CreateVideo
                    fas.Position=pos;
                end
                hold off
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,dhdt,CtrlVar);
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                
                title(sprintf('dhdt at t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                title(colorbar,'(m)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
                caxis([-10 10])
            end
            
            
            if ~isempty(strfind(plots,'-h-'))
                if ~exist('fh','var') || ~ishandle(fh)
                    fh=figure;
                else
                    figure(fh)
                end
                
                if CreateVideo
                    fh.Position=pos;
                    %hfig.Position=pos;
                end
                hold off
                
                %fh.NextPlot='replacechildren';  % slows things down
                
                
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,h,CtrlVar);
                hold on ;
                %plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                
                if ~isempty(strfind(CtrlVar.Experiment,'Jenkins'))
                    title(sprintf('ice thickness at t=%-g. Tw=%-g (yr)',CtrlVar.time,Tw)) ;
                else
                    title(sprintf('ice thickness at t=%-g (yr)',CtrlVar.time)) ;
                end
                
                
                xlabel('xps (km)') ; ylabel('yps (km)')
                title(colorbar,'(m)')
                %                 if PlotMinThickLocations
                %                     plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                %                 end
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
                
                caxis([100 800])
            end
            
            
            if ~isempty(strfind(plots,'-sbB-'))
                if ~exist('fsbB','var') || ~ishandle(fsbB)
                    fsbB=figure;
                else
                    figure(fsbB)
                end
                if CreateVideo
                    fsbB.Position=pos;
                end
                hold off
                
                AspectRatio=1; ViewAndLight=[-45 35 -45 50];
                [TRI,DT]=Plot_sbB(CtrlVar,MUA,s,b,B,TRI,DT,AspectRatio,ViewAndLight);
                
                title(sprintf('t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
                
            end
            
            if ~isempty(strfind(plots,'-MeltNodes-'))
                if ~exist('MN','var') || ~ishandle(MN)
                    MN=figure;
                else
                    figure(MN)
                end
                if CreateVideo
                    MN.Position=pos;
                end
                hold off
                
                CtrlVar.MeltNodesDefinition='Element-Wise';
                [MeltNodes,NotMeltNodes]=SpecifyMeltNodes(CtrlVar,MUA,GF);
                PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
                hold on
                plot(x(MeltNodes)/CtrlVar.PlotXYscale,y(MeltNodes)/CtrlVar.PlotXYscale,'or',...
                    'MarkerSize',3,'MarkerFaceColor','r')
                hold on
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                xlabel(CtrlVar.PlotsXaxisLabel) ; ylabel(CtrlVar.PlotsYaxisLabel) ;
                title(['Melt nodes: time=',num2str(CtrlVar.time)])
                
                
                title(sprintf('Melt Nodes at t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
            end
            
            
            if ~isempty(PlotRegion)
                SetRegionalPlotAxis(PlotRegion);
            else
                if ~isnan(PlotArea)
                    axis(PlotArea)
                end
            end
           
            %drawnow
            
            if CreateVideo
                iFrame=iFrame+1;
                Frame = getframe(gcf);
                %Frame = hardcopy(hFig, '-opengl', '-r0');
                writeVideo(vidObj,Frame);
                hold off
            else
                
                if ~strcmpi(usrstr,'C')
                    prompt = 'For next plot press RET. To go back enter B. To exit enter E: ';
                    usrstr = input(prompt,'s');
                    if strcmpi(usrstr,'E')
                        break
                    elseif strcmpi(usrstr,'B')
                        iFile=IIlist(k-iback)-1;
                        iback=iback+1;
                    else
                        iback=0;
                        k=k+1; IIlist(k)=iFile;
                    end
                end
            end
            PlotArea=axis;
        end
        
%{        
        reply = input('Do you want more? Y/N [Y]:','s');
        if isempty(reply)
            reply = 'Y';
        else
            break
        end
%}        
        iFile=iFile+1;
        
    end
    
    if CreateVideo
        close(vidObj);
        fprintf('\n video file closed \n')
    end
  
    close all
end

cd(CurDir)



