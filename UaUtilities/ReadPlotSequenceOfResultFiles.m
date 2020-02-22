function DataCollect=ReadPlotSequenceOfResultFiles(FileNameSubstring,PlotType,PlotScreenPosition,N,AxisLimits)


%%
% A simple plotting routine to plot a sequence of results files.
% This is just a template and will needed to be adjusted to suite a given situation
%%
% Approach:
% 1) look for files containing a given substring in a ./ResultsFiles subdirectory
% 2) assume the first N letters of the filenames contain the time (multiplied by 100)
% 3) loop over results files, read the results and create plots
%
%
% Example:
%
%   DataCollect=ReadPlotSequenceOfResultFiles("-Ice1r-","-collect-");
%   figure ; plot(DataCollect.time, DataCollect.VAF/1e9,'-or'); xlabel("time (yr)") ; ylabel(" VAF (Gt)")
%   figure ; plot(DataCollect.time, DataCollect.GroundedArea/1e6,'-or'); xlabel("time (yr)") ; ylabel(" Grounded area (km^2)")
%
% 
% Example:
%
%    Ice1r=ReadPlotSequenceOfResultFiles("-Ice1r-","-collect-");
%    Ice1ra=ReadPlotSequenceOfResultFiles("-Ice1ra-","-collect-");
%    figure ; plot(Ice1r.time,Ice1r.GroundedArea/1e6,'-o');
%    hold on ; plot(Ice1ra.time,Ice1ra.GroundedArea/1e6,'-x');
%    legend("Ice1r","ice1ra")
%    xlabel("time (yr)") ; ylabel(" VAF (Gt)")
%
% Example:
%
%    Ice1r=ReadPlotSequenceOfResultFiles("-Ice1r-","-mesh-");
%
%
%    Ice2r=ReadPlotSequenceOfResultFiles("-Ice2r-Tsai-","-mesh-speed-s-ab-");
%
%% Parameters

narginchk(0,7)
nargoutchk(0,1)


if nargin==0 || isempty(FileNameSubstring)
    FileNameSubstring="";
end

if nargin<2 || isempty(PlotType)
    
    PlotType="-mesh-";  %  specify the type of plot to create, see below, modify and expand as needed
    PlotType="-h-";
    PlotType="-sbB-";
    PlotType="-MeltNodes-";
    PlotType="-ab-";
    PlotType="-dhdt-";
    PlotType="-ubvb-";
    PlotType="-log10(BasalSpeed)-";
    %PlotType="-VAF-";
    PlotType="-collect-";
    PlotType="-mesh-speed-s-ab-";
end


if nargin<3 || isempty(PlotScreenPosition)
    PlotScreenPosition=[40 40 2300 1800]; % positon of figure on screen, adjust this to your own screen resolution
    PlotScreenPosition=[40 40 2000 1600]; %
end

if nargin<4 || isempty(N)
    N=7;
end

PlotTimeInterval=1;                     % model time interval between creation of plots
PlotTimeMax=1e10;

if contains(PlotType,"-collect-")
    CreateVideo=0;
else
    CreateVideo=1;
end

PlotRegion=[];
PlotMinThickLocations=true;

%%
CurDir=pwd;

if CreateVideo
    
    vidObj = VideoWriter("VideoResultsFile"+FileNameSubstring+PlotType);
    vidObj.FrameRate=1;   % frames per sec
    open(vidObj);
end

% assume results files have been saved in a subdirectory named 'ResultsFiles' (modify as needed)
% cd ./ResultsFiles/
SearchString="*"+FileNameSubstring+"*.mat";
SearchString=replace(SearchString,"**","*");
list=dir(SearchString);   % get names of all .mat files containing a given substring



% cd ..   ; % go back up into working directory


CtrlVar.QuiverSameVelocityScalingsAsBefore=false;
nFiles=length(list);
iFile=1; iFrame=1;
iCount=0;
DataCollect=[];
DataCollect.FileNameSubstring=FileNameSubstring;

while iFile<=nFiles   % loop over files
    
    
    time=str2double(list(iFile).name(1:N))/100;  % get the model time, assuming that the first N letters of filename are the model time*100
    %time=str2double(list(iFile).name(1:4));
    if mod(time,PlotTimeInterval)==0 && time<=PlotTimeMax   % only do plots at given time intervals and up to a max time specifed
        
        try   % go back into subdirectory containing result files and load one result file
          
            load(list(iFile).name)
           
            fprintf(' %s \n ',list(iFile).name)
        catch
            fprintf('could not load %s \n ',list(iFile).name)
        end
        
        CtrlVar.PlotXYscale=1000;
        [GLgeo,GLnodes,GLele]=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar);
        
        TRI=[]; DT=[]; xGL=[];yGL=[];
        x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2);
        ih=F.h<=CtrlVar.ThickMin;
        
        iCount=iCount+1;
        DataCollect.time(iCount)=CtrlVar.time;
        
        % if creaing 4 subplots, try to arrange them based on the aspect ratio
        XYratio=(max(MUA.coordinates(:,2))-min(MUA.coordinates(:,2)))/(max(MUA.coordinates(:,1))-min(MUA.coordinates(:,1)));
        
        if XYratio>0.5 || XYration<1.5
            nPx=2 ; nPy=2;
        elseif XYratio<=0.5
            nPx=1; nPy=4;
        else
            nPx=4; nPy=1;
        end
        
        if nargi<
        
        switch PlotType
            
            case "-collect-"
                
                % collect data across all result files
                
                [VAF,IceVolume,GroundedArea]=CalcVAF(CtrlVar,MUA,F.h,F.B,F.S,F.rho,F.rhow,F.GF);
                
                DataCollect.VAF(iCount)=VAF.Total ;
                DataCollect.GroundedArea(iCount)=GroundedArea.Total;
                DataCollect.IceVolume(iCount)=IceVolume;
                
                
                %%
                
            case "-mesh-speed-s-ab-"
                
                f4=FindOrCreateFigure("-mesh-speed-s-ab-",PlotScreenPosition);
                clf(f4)
                hold off
                
                subplot(nPx,nPy,1)
                hold off
                PlotMuaMesh(CtrlVar,MUA);
                title(sprintf('t=%-g (yr)  #Ele=%-i, #Nodes=%-i, #nod=%-i',time,MUA.Nele,MUA.Nnodes,MUA.nod))
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r','LineWidth',2);
                xlabel('x (km)') ; ylabel('y (km)') ;
                axis equal tight 
                hold off
                
                subplot(nPx,nPy,2)
                speed=sqrt(F.ub.*F.ub+F.vb.*F.vb);
                hold off
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,speed); title(sprintf('speed at t=%g',time))
                %QuiverColorGHG(MUA.coordinates(:,1),MUA.coordinates(:,2),ub,vb,CtrlVar);
                hold on
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r','LineWidth',2);
                xlabel('x (km)') ; ylabel('y (km)') ; title(cbar,'(m/yr)')
                hold off
                
                subplot(nPx,nPy,3)
                hold off
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.s);   title(sprintf('surface at t=%g',time))
                hold on
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r','LineWidth',2);
                xlabel('x (km)') ; ylabel('y (km)') ; title(cbar,'(m/yr)')
                
                subplot(nPx,nPy,4)
                hold off
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.ab);   title(sprintf('Basal melt at t=%g',time))
                hold on
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r','LineWidth',2);
                xlabel('x (km)') ; ylabel('y (km)') ; title(cbar,'(m/yr)')
                hold off
                
                
                sgtitle(FileNameSubstring)
                
                %%
                
                
            case '-mesh-'
                
                fmesh=FindOrCreateFigure('Mesh',PlotScreenPosition);
                
                PlotMuaMesh(CtrlVar,MUA);
                title(sprintf('t=%-g (yr)  #Ele=%-i, #Nodes=%-i, #nod=%-i',time,MUA.Nele,MUA.Nnodes,MUA.nod))
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
                hold off
                
                
                
                %%
                
            case '-ubvb-'
                % plotting horizontal velocities
                %%
                
                fubvb=FindOrCreateFigure('fubvb',PlotScreenPosition);
                
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
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'k');
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
                
                
                SurfSpeed=sqrt(F.ub.*F.ub+F.vb.*F.vb);
                
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(SurfSpeed),CtrlVar);
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'k');
                
                title(sprintf('log_{10}(Basal speed) t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                title(colorbar,'log_{10}(m/yr)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                caxis([1 4])
                %%
                
            case "-VAF-"
                
                
                figlogSpeed=FindOrCreateFigure('VAF',PlotScreenPosition);
                
                [VAF,IceVolume,GroundedArea]=CalcVAF(CtrlVar,MUA,F.h,F.B,F.S,F.rho,F.rhow,F.GF);
                
                if iCount==0
                    vaf=[];
                end
                
                
                vaf.value(iCount)=VAF.Total ;
                vaf.time(iCount)=CtrlVar.time;
                hold off
                plot(vaf.time,vaf.value/1e9,'or')
                xlabel(' time (yr)') ; ylabel(' VAF (Gt)')
                
                
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
                
                F.ab(F.GF.node>0.5)=NaN;
                
                hold off
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,F.ab,CtrlVar);
                
                hold on
                
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'k');
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
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'k');
                
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
                
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'k');
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
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],'LineWidth',2);
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
                
                CtrlVar.MeltNodesDefinition='Node-Wise';
                [MeltNodes,NotMeltNodes]=SpecifyMeltNodes(CtrlVar,MUA,F.GF,GLgeo,GLnodes,GLele);
                PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
                hold on
                plot(x(MeltNodes)/CtrlVar.PlotXYscale,y(MeltNodes)/CtrlVar.PlotXYscale,'o',...
                    'MarkerSize',3,'MarkerFaceColor','g')
                hold on
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
                xlabel(CtrlVar.PlotsXaxisLabel) ; ylabel(CtrlVar.PlotsYaxisLabel) ;
                title(['Melt nodes: time=',num2str(CtrlVar.time)])
                
                
                title(sprintf('Melt Nodes at t=%-g (yr)',time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                
                %                 if PlotMinThickLocations
                %                     plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                %                 end
                
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
        
        if ~contains(PlotType,"-collect-")
            PlotArea=axis;
        end
    end
    iFile=iFile+1;
    
end

if CreateVideo
    close(vidObj);
    fprintf('\n video file closed \n')
end

close all

cd(CurDir)

end

