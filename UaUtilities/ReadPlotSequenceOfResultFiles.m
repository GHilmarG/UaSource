%%
% A simple plotting routine to plot a sequence of results files.
% This is just a template and will needed to be adjusted to suite a given situation
%%
% Approach:
% 1) look for files containing a given substring in a ./ResultsFiles subdirectory
% 2) assume the first 10 letters of the filenames contain the time (multiplied by 100)
% 3) loop over results files, read the results and create plots
%

%% Parameters
FileNameSubstring='UserOutputFile-prognostic-n3-m5-MatlabOptimization-Nod3-I-Adjoint-Cga1-Cgs1-Aga1-Ags1-0-0-logAGlenlogC';
PlotTimeInterval=1;                     % model time interval between creation of plots
PlotTimeMax=1e10;
PlotType='-mesh-';  %  specify the type of plot to create, see below, modify and expand as needed
PlotType='-log10(BasalSpeed)-';
PlotType='-dhdt-';
PlotType='-ubvb-';
PlotType='-ab-';
PlotType='-h-';
PlotType='-sbB-';
PlotScreenPosition=[40 40 2300 1800]; % positon of figure on screen, adjust this to your own screen resolution
PlotScreenPosition=[40 40 2000 1600]; %
%pos=[200 50 1200 900];
CreateVideo=1;
PlotRegion=[];
PlotMinThickLocations=true;

%%
CurDir=pwd;

if CreateVideo
    
    vidObj = VideoWriter(['VideoResultsFile',PlotType,'.mp4']);
    vidObj.FrameRate=2;
    %vidObj.FrameRate=1;   % frames per sec
    open(vidObj);
end

% assume results files have been saved in a subdirectory named 'ResultsFiles' (modify as needed)
cd ./ResultsFiles/
list=dir(['*',FileNameSubstring,'*.mat']);   % get names of all .mat files containing a given substring
cd ..   ; % go back up into working directory


CtrlVar.QuiverSameVelocityScalingsAsBefore=false;
nFiles=length(list);
iFile=1; iFrame=1;

while iFile<=nFiles   % loop over files
    
    
    time=str2double(list(iFile).name(1:10))/100;  % get the model time, assuming that the first 10 letters of filename are the model time*100
    %time=str2double(list(iFile).name(1:4));
    if mod(time,PlotTimeInterval)==0 && time<=PlotTimeMax   % only do plots at given time intervals and up to a max time specifed
        
        try   % go back into subdirectory containing result files and load one result file
            cd ./ResultsFiles/
            load(list(iFile).name)
            cd ..
            fprintf(' %s \n ',list(iFile).name)
        catch
            fprintf('could not load %s \n ',list(iFile).name)
        end
        
        CtrlVar.PlotXYscale=1000;
        GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
        TRI=[]; DT=[]; xGL=[];yGL=[];
        x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2);
        ih=F.h<=CtrlVar.ThickMin;
        
        
        switch PlotType
            
            case '-mesh-'
                
                
                if ~exist('fmesh','var') || ~ishandle(fmesh)
                    fmesh=figure;
                else
                    figure(fmesh)
                end
                
                if CreateVideo
                    fmesh.Position=PlotScreenPosition;
                end
                
                
                hold off
                PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
                title(sprintf('t=%-g (yr)  #Ele=%-i, #Nodes=%-i, #nod=%-i',time,MUA.Nele,MUA.Nnodes,MUA.nod))
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
                
                
                
            case '-ubvb-'
                % plotting horizontal velocities
                %%
                if ~exist('fubvb','var') || ~ishandle(fubvb)
                    fubvb=figure;
                else
                    figure(fubvb)
                end
                
                if CreateVideo
                    fubvb.Position=PlotScreenPosition;
                end
                
                
                hold off
                N=1;
                %speed=sqrt(ub.*ub+vb.*vb);
                %CtrlVar.MinSpeedWhenPlottingVelArrows=0; CtrlVar.MaxPlottedSpeed=max(speed); %
                CtrlVar.VelPlotIntervalSpacing='log10';
                %CtrlVar.VelColorMap='hot';
                %CtrlVar.RelativeVelArrowSize=10;
                CtrlVar.QuiverColorSpeedLimits=[1 5000];
                
                PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,'b')
                hold on
                QuiverColorGHG(x(1:N:end),y(1:N:end),F.ub(1:N:end),F.vb(1:N:end),CtrlVar);
                CtrlVar.QuiverSameVelocityScalingsAsBefore=true;
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                title(sprintf('(ub,vb) t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                %%
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                end
                
                
                
                
                
            case '-log10(BasalSpeed)-'
                %%
                %us=ub+ud;  vs=vb+vd;
                if ~exist('fab','var') || ~ishandle(fab)
                    fab=figure;
                else
                    figure(fab)
                end
                if CreateVideo
                    fab.Position=PlotScreenPosition;
                end
                SurfSpeed=sqrt(F.ub.*F.ub+F.vb.*F.vb);
                hold off
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(SurfSpeed),CtrlVar);
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                
                title(sprintf('log_{10}(Basal speed) t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                title(colorbar,'log_{10}(m/yr)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                caxis([1 4])
                %%
                
                
            case '-ab-'
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
                    fab.Position=PlotScreenPosition;
                    %fab.NextPlot='replacechildren';  % slows things down
                end
                
                %ax.NextPlot='replace';
                
                F.ab(GF.node>0.5)=NaN;
                
                hold off
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,F.ab,CtrlVar);
                
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
                end
                caxis([min(F.ab) max(F.ab)]) ;
                
                
                
            case '-dhdt-'
                if ~exist('fas','var') || ~ishandle(fas)
                    fas=figure;
                else
                    figure(fas)
                end
                if CreateVideo
                    fas.Position=PlotScreenPosition;
                end
                hold off
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,F.dhdt,CtrlVar);
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
                caxis([-15 15])
                
                
                
            case '-h-'
                
                if ~exist('fh','var') || ~ishandle(fh)
                    fh=figure;
                else
                    figure(fh)
                end
                
                if CreateVideo
                    fh.Position=PlotScreenPosition;
                end
                hold off
                
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,F.h,CtrlVar);
                hold on ;
                
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
                title(sprintf('ice thickness at t=%-g (yr)',CtrlVar.time)) ;
                title(colorbar,'(m)')
                
                xlabel('xps (km)') ; ylabel('yps (km)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                end
                
            case '-sbB-'
                if ~exist('fsbB','var') || ~ishandle(fsbB)
                    fsbB=figure;
                else
                    figure(fsbB)
                end
                if CreateVideo
                    fsbB.Position=[100 100 2000 1300];
                end
                hold off
                
                CtrlVar.PlotGLs=0;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'LineWidth',2);
                Fs=scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.s,'natural');
                zGL=Fs(xGL,yGL);
                speed=sqrt(F.ub.*F.ub+F.vb.*F.vb);
                minSpeed=0 ; maxSpeed=4000;
                speed(speed<(minSpeed+eps))=minSpeed;
                speed(speed>maxSpeed)=maxSpeed;
                
                sCol=jet(numel(F.s));
                sCol=parula(numel(F.s));
                ColorIndex=Variable2ColorIndex(speed,[minSpeed maxSpeed]);
                sCol(:,:)=sCol(ColorIndex,:);
                sCol(speed<minSpeed,:)=sCol(speed<minSpeed,:)*0+[0.9 0.9 0.9];
                %bCol=jet(numel(s))*0+[0 0 1];
                
                brown=[0.59 0.29 0.00];
                BCol=jet(numel(F.s))*0+brown; % [0.9 0.9 0.9];
                bCol=BCol;
                AspectRatio=1; ViewAndLight=[-45 25 -45 50];
                LightHandle=[];
                [TRI,DT]=Plot_sbB(CtrlVar,MUA,F.s,F.b,F.B,TRI,DT,AspectRatio,ViewAndLight,LightHandle,sCol,bCol,BCol);
                xlabel('x (km)' ) ; ylabel('y (km)' )
                hold on
                plot3(xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,zGL,'color','r','LineWidth',2)
                %colormap(jet)
                colormap(parula)
                
                cbar=colorbar;
                
                for I=1:numel(cbar.TickLabels)
                    cbar.TickLabels{I}=str2double(cbar.TickLabels{I})*(max(speed)-min(speed))+min(speed);
                end
                
                title(cbar,'(m/yr)')
                
                fsbB.CurrentAxes.DataAspectRatio=[1 1 20];
                view(-55,35)
                %axis([-1700 -1300 -600 -100 -2000  3000]) % small
                
                camproj('perspective')
                title(sprintf('t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                
                
            case '-MeltNodes-'
                
                if ~exist('MN','var') || ~ishandle(MN)
                    MN=figure;
                else
                    figure(MN)
                end
                if CreateVideo
                    MN.Position=PlotScreenPosition;
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
                end
                
                
            case '-Bab-'
                
                if ~isempty(PlotRegion)
                    SetRegionalPlotAxis(PlotRegion);
                else
                    if ~isnan(PlotArea)
                        axis(PlotArea)
                    end
                end
                
                %drawnow
                Sub1.Position=[0.05 0.55 0.35 0.43];
                Sub2.Position=[0.55 0.55 0.35 0.43];
                Sub3.Position=[0.05 0.05 0.35 0.43];
                Sub4.Position=[0.55 0.05 0.35 0.43];
                
        end
        
        
        if CreateVideo
            iFrame=iFrame+1;
            Frame = getframe(gcf);
            %Frame = hardcopy(hFig, '-opengl', '-r0');
            writeVideo(vidObj,Frame);
            hold off
            
        end
        PlotArea=axis;
        
    end
    iFile=iFile+1;
    
end

if CreateVideo
    close(vidObj);
    fprintf('\n video file closed \n')
end

close all

cd(CurDir)



