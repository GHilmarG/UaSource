function DataCollect=ReadPlotSequenceOfResultFiles(varargin)


%%
% A simple plotting routine to plot a sequence of results files. This is just a template
% and will needed to be adjusted to suite a given situation.
%
%%
% Approach:
%
% 1) look for files containing a given substring in a current directory
%
% 2) assume the first N letters of the filenames contain the time (multiplied by 100). By defaul N=7.
%
% 3) loop over results files, read the results and create plots
%
% Example:
%
%
% Read in all .mat file in current directory and use default plot settings for plotting
%
%   ReadPlotSequenceOfResultFiles();
%
% Example:
%
% Read only those .mat files that have FileNameSubstring as a part of their name, collect
% some data and plot VAF and grounded area as function of time.
%
%   DataCollect=ReadPlotSequenceOfResultFiles(FileNameSubstring="-Ice1r-",PlotType="-collect-");
%   figure ; plot(DataCollect.time, DataCollect.VAF/1e9,'-or'); xlabel("time (yr)") ;
%   ylabel(" VAF (Gt)") figure ; plot(DataCollect.time,DataCollect.GroundedArea/1e6,'-or');
%   xlabel("time (yr)") ; ylabel(" Grounded area(km^2)")
%
%
% Examples:
%
%
%   ReadPlotSequenceOfResultFiles(FileNameSubstring="ResultsFile",PlotTimestep=10)
%
% Read only those .mat files that have FileNameSubstring as a part of their name, set the
% axis limits, and plot a file every 0.1 time units.
%
%    ReadPlotSequenceOfResultFiles("FileNameSubstring","Forward","AxisLimits",[480 494 -2296 -2280],"PlotTimestep",10) ;
%
%
%
%   ReadPlotSequenceOfResultFiles("FileNameSubstring","TravellingFront-1dAnalyticalIceShelf","PlotType","-1dIceShelf-","PlotTimestep",1,"PlotTimeInterval",[0 200])
%   ReadPlotSequenceOfResultFiles("FileNameSubstring","TravellingFront-1dAnalyticalIceShelf-MBice0-SUPGtaus-Adapt1","PlotType","-1dIceShelf-","PlotTimestep",1,"PlotTimeInterval",[0 200])
%   ReadPlotSequenceOfResultFiles("FileNameSubstring","Ex-Calving-1dIceShelf-MBice0-SUPGtaus-Adapt1","PlotType","-1dIceShelf-","PlotTimestep",1,"PlotTimeInterval",[0 200])
%
%    ReadPlotSequenceOfResultFiles("FileNameSubstring","Ex-Calving-1dIceShelf-MBice0-SUPGtaus-Adapt1","PlotType","-Level Set-","PlotTimestep",10,"PlotTimeInterval",[0 2000],"PlotScreenPosition",[100 650 1100 570])
%
%    ReadPlotSequenceOfResultFiles("FileNameSubstring","Forward-Transient-Nod3-rCW-N0-M-Cga1-Cgs10000-Aga1-Ags10000-logAGlenlogC-uv-","PlotType","-Level Set-","PlotTimestep",1)
%
%    ReadPlotSequenceOfResultFiles("FileNameSubstring","LevelSetWithMeltFeedback","PlotType","-1dIceShelf-","PlotTimestep",1)
%   
%
%% Parse inputs


defaultPlotType="-mesh-speed-calving-level set-";
allowedPlotTypes = {'-mesh-speed-s-ab-','-mesh-speed-calving-level set-','-mesh-',...
    '-h-','-sbB-','-dhdt-','-log10(BasalSpeed)-','-VAF-',...
    '-1dIceShelf-','-hPositive-','-Level Set-',...
    '-collect-'};


IP = inputParser;

addParameter(IP,"FileNameSubstring","",@isstring);
addParameter(IP,"VideoFileName","VideoFile",@isstring);
addParameter(IP,"AxisLimits",NaN,@isnumeric);
addParameter(IP,"N",7,@isnumeric);
addParameter(IP,"PlotType",defaultPlotType,@(x) any(validatestring(x,allowedPlotTypes)));
addParameter(IP,"PlotScreenPosition",NaN,@isnumeric);
addParameter(IP,"PlotTimeInterval",[0 1e10],@isnumeric);
addParameter(IP,"PlotTimestep",1,@isnumeric);

addParameter(IP,"DataToBeCollected","",@isstring);


parse(IP,varargin{:});

FileNameSubstring=IP.Results.FileNameSubstring ;
VideoFileName=IP.Results.VideoFileName+FileNameSubstring ;
PlotType=IP.Results.PlotType ;
N=IP.Results.N ;
AxisLimits=IP.Results.AxisLimits;
PlotScreenPosition=IP.Results.PlotScreenPosition;
PlotTimeInterval=IP.Results.PlotTimeInterval;
PlotTimestep=IP.Results.PlotTimestep;
DataToBeCollected=IP.Results.DataToBeCollected; 
%%

if PlotType=="-collect-" && DataToBeCollected==""  %
    DataToBeCollected="-VAF-Grounded Area-" ;
end

if  DataToBeCollected~=""
    PlotType="-collect-" ;  % make sure that if DataToBeCollected is not empty, the PlotType is collect.
end

if isnan(PlotScreenPosition)
    
    ScreenSize=get(0,'ScreenSize')+10;
    HW=min(ScreenSize(3),ScreenSize(4));
    PlotScreenPosition=[5 5 HW-10 HW-10];  % creates a square plot window
end


if contains(PlotType,"-collect-")
    CreateVideo=0;
else
    CreateVideo=1;
end

PlotRegion=[];
PlotArea=[];
PlotMinThickLocations=true;

%%
CurDir=pwd;

if CreateVideo
    
 
    VideoFileName = replace(VideoFileName,".mat","");
    VideoFileName = replace(VideoFileName," ","");
    vidObj = VideoWriter(VideoFileName);
    vidObj.FrameRate=10;   % frames per sec
    open(vidObj);
end

% assume results files have been saved in a subdirectory named 'ResultsFiles' (modify as needed)
% cd ./ResultsFiles/
if ~contains(FileNameSubstring,".mat")
    SearchString="*"+FileNameSubstring+"*.mat";
else
    SearchString="*"+FileNameSubstring;
end
SearchString=replace(SearchString,"**","*");
list=dir(SearchString);   % get names of all .mat files containing a given substring



% cd ..   ; % go back up into working directory


CtrlVar.QuiverSameVelocityScalingsAsBefore=false;
nFiles=length(list);
iFile=1; iFrame=1;
iCount=0;
DataCollect=[];
DataCollect.FileNameSubstring=FileNameSubstring;
LSFScale=[];

while iFile<=nFiles   % loop over files
    
    
    time=str2double(list(iFile).name(1:N))/100;  % get the model time, assuming that the first N letters of filename are the model time*100
    %time=str2double(list(iFile).name(1:4));
    if mod(time,PlotTimestep)==0 && time<=PlotTimeInterval(2) && time>=PlotTimeInterval(1)   % only do plots at given time intervals and up to a max time specifed
        
        try   % go back into subdirectory containing result files and load one result file
            
            load(list(iFile).name)
            
            fprintf(' %s \n ',list(iFile).name)
            
            % Just in case the results files are old an from the time when F.GF and F.x and F.y were not defined as F fields. 
            if ~isfield(F,"GF") && exist("GF","var")
                F.GF=GF;
            end
            if  ~isfield(F,"x")
                F.x=MUA.coordinates(:,1);
                F.y=MUA.coordinates(:,2);
            end
            
        catch
            fprintf('could not load %s \n ',list(iFile).name)
        end
        
        CtrlVar.PlotXYscale=1000;
        [GLgeo,GLnodes,GLele]=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar);
        
        TRI=[]; DT=[]; xGL=[];yGL=[];
        x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2);
        ih=F.h<=CtrlVar.ThickMin;
        
        
        
        
        % if creaing 4 subplots, try to arrange them based on the aspect ratio
        
        xmin=min(MUA.coordinates(:,1)) ;
        xmax=max(MUA.coordinates(:,1)) ;
        ymin=min(MUA.coordinates(:,2)) ;
        ymax=max(MUA.coordinates(:,2)) ;
        
        XYratio=(xmax-xmin)/(ymax-ymin) ;
        
        if XYratio>0.5 && XYratio<1.5
            nPx=2 ; nPy=2;
        elseif XYratio<=0.5
            nPx=1; nPy=4;
        else
            nPx=4; nPy=1;
        end
        
       
        if ~isnan(AxisLimits) ; axis(AxisLimits) ; end
        
        switch PlotType
            
            case "-collect-"
                
                % collect data across all result files
                
                if ~isfield(DataCollect,'time')
                    DataCollect.time=zeros(nFiles,1)+NaN;
                    DataCollect.VAF=zeros(nFiles,1)+NaN;
                    DataCollect.GroundedArea=zeros(nFiles,1)+NaN;
                    DataCollect.IceVolume=zeros(nFiles,1)+NaN;
                    
                    DataCollect.IceVolume=zeros(nFiles,1)+NaN;
                    
                    DataCollect.xcMax=zeros(nFiles,1)+NaN;
                    DataCollect.xcMin=zeros(nFiles,1)+NaN;
                    DataCollect.xcMean=zeros(nFiles,1)+NaN;
                    DataCollect.Lx=zeros(nFiles,1)+NaN;
                    DataCollect.xc=zeros(nFiles,1)+NaN;
                    
                    if contains(DataToBeCollected,"-GL Flux-")
                        DataCollect.qGL=zeros(nFiles,1)+NaN;
                    end
                    
                end

                
                
                [VAF,IceVolume,GroundedArea]=CalcVAF(CtrlVar,MUA,F.h,F.B,F.S,F.rho,F.rhow,F.GF);
                
                iCount=iCount+1;
                DataCollect.time(iCount)=CtrlVar.time;
                DataCollect.VAF(iCount)=VAF.Total ;
                DataCollect.GroundedArea(iCount)=GroundedArea.Total;
                DataCollect.IceVolume(iCount)=IceVolume.Total;
                
                DataCollect.Lx(iCount)=max(F.x(F.h>10)) ;
                
                if contains(DataToBeCollected,"-GL Flux-")
                    [qGL,qGLx,qGLy,~,~,~,~,~,GLgeo]=FluxAcrossGroundingLine(CtrlVar,MUA,F.GF,F.ub,F.vb,F.ud,F.vd,F.h,F.rho) ;
                    DataCollect.qGL(iCount)=sum(qGL);
                end
                 
                 
                if isfield(F,"LSF") &&  ~isempty(F.LSF)
                    [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F.LSF,0);
                    
                    DataCollect.xcMax(iCount)=max(xc,[],'omitnan') ;
                    DataCollect.xcMin(iCount)=min(xc,[],'omitnan') ;
                    DataCollect.xcMean(iCount)=mean(xc,'omitnan') ;
                end
                
                %%
                
            case "-hPositive-"
                
                figV=FindOrCreateFigure('-hPositive-',[100 100 1500 1500]) ;
                [TotalIceVolume,ElementIceVolume]=CalcIceVolume(CtrlVar,MUA,F.h);
                I=F.h<=(CtrlVar.ThickMin+100*eps);
                Reactions=CalculateReactions(CtrlVar,MUA,BCs,l);
                
                subplot(2,2,1)
                TRI = delaunay(x,y);
                trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,F.h,100*F.s+100,'EdgeColor','none') ; hold on
                
                plot3(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,F.h(I),'ko');
                
                view(45,25); lightangle(-45,10) ; lighting phong ;
                xlabel('y (km)') ; ylabel('x (km)') ; zlabel('ice thickness (m)') ;
                %colorbar ; title(colorbar,'(m)')
                colormap(flipud(othercolor('RdYlBu_11b',2000)))
                hold on
                
                nThickConstraints=numel(find(I));
                title(sprintf(' Volume=%#8.3f (m^3)  #Contraints=%i',TotalIceVolume/1e9,nThickConstraints))
                axis normal
                zlim([0 110])
                %axis equal ; tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)*CtrlVar.PlotXYscale/80]); axis tight
                hold off
                
                subplot(2,2,2)
                hold off
                %                     if ~isempty(Reactions.h)
                %                         [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,Reactions.h./F.rho) ;
                %                         title(cbar,'Reactions (m)')
                %                         colormap(flipud(othercolor('RdYlBu_11b',2000))) ;
                %                     end
                
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.h) ; title(cbar,'Ice thickness (m)')
                caxis([-0.1 1])
                colormap(flipud(othercolor('RdYlBu_11b',2000))) ;
                
                if ~isnan(AxisLimits) ; axis(AxisLimits) ; end
                xlabel('x (km)') ; ylabel('y (km)')
                hold on
                plot(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,'wo');
                hold off
                
                subplot(2,2,3)
                hold on
                yyaxis left
                
                plot(CtrlVar.time,TotalIceVolume/1e9,'o') ;
                ylabel('Ice Volume (m^3)')
                ylim([3.7e3 4e3])
                
                rDist=(MUA.coordinates(:,1).^2+MUA.coordinates(:,2).^2);
                [temp,I]=min(rDist);
                hCentre=F.s(I)-F.b(I);
                yyaxis right
                plot(CtrlVar.time,hCentre,'o') ;
                ylabel("Centre ice thickness (m)")
                ylim([0 10]) ; % ylim([0 101])
                
                
                xlabel('time (yr)') ;
                xlim(PlotTimeInterval)
                
                subplot(2,2,4)
                hold off ;
                
                if isempty(Reactions.h)
                    Reactions.h=zeros(MUA.Nnodes,1);
                end
                
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.h.*Reactions.h./F.rho) ;
                caxis([min(F.h.*Reactions.h./F.rho)-eps max(F.h.*Reactions.h./F.rho)+eps])
                
                colormap(flipud(othercolor('RdYlBu_11b',1000))) ;
                ModifyColormap ;
                title(cbar,'$\lambda h$ ($m^2$)','Interpreter','latex')
                
                
                hold on
                
                plot(x(BCs.hPosNode)/CtrlVar.PlotXYscale,y(BCs.hPosNode)/CtrlVar.PlotXYscale,'*b');
                if numel(BCs.hPosNode) >0
                    scale=1000;
                    
                    Ir=Reactions.h>0; scatter(x(Ir)/CtrlVar.PlotXYscale,y(Ir)/CtrlVar.PlotXYscale,+scale*Reactions.h(Ir)./F.rho(Ir),'r')
                    Ir=Reactions.h<0; scatter(x(Ir)/CtrlVar.PlotXYscale,y(Ir)/CtrlVar.PlotXYscale,-scale*Reactions.h(Ir)./F.rho(Ir),'b')
                    
                end
                PlotMuaMesh(CtrlVar,MUA);
                
                xlabel('x (km)') ; ylabel('y (km)')
                
                if ~isnan(AxisLimits) ; axis(AxisLimits) ; end
                sgtitle(sprintf("Cubic element.  time=%3.0f",CtrlVar.time))
                
                
            case "-Level Set-"
                
                if ~isempty(F.LSF)
                    
                    
                    fMeshLSF=FindOrCreateFigure("Mesh and LSF");
                    
                    
                    % fMeshLSF.Position=[100 650 1100 570] ;
                    
                    %[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.h) ; title(cbar,'Thickness (m)') ; caxis([0 50])
                    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.LSF/1000) ; title(cbar,'LSF (km)') ; caxis([-10 10])
                    colormap(flipud(othercolor('RdYlBu_11b',2000)))
                    
                    hold on
                    PlotMuaMesh(CtrlVar,MUA,'k'); hold on
                    
                    
                    if ~isempty(F.LSFMask)
                        plot(MUA.coordinates(F.LSFMask.NodesOut,1)/1000,MUA.coordinates(F.LSFMask.NodesOut,2)/1000,'*r')
                    end
                    
                    
                    [xGL,yGL]=PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],'r','LineWidth',2);
                    if ~isempty(xGL)
                        Temp=fMeshLSF.CurrentAxes.Title.String;
                        fMeshLSF.CurrentAxes.Title.String=[Temp(:)',{"Grounding line in red"}];
                    end
                    
                    if ~isempty(F.LSF) && CtrlVar.LevelSetMethod
                        hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,'w','LineWidth',2) ;
                        Temp=fMeshLSF.CurrentAxes.Title.String;
                        fMeshLSF.CurrentAxes.Title.String=[Temp(:)',{"Level-set zero line in white"}];
                        
                        [xf,yf]=CalcMuaFieldsContourLine(CtrlVar,MUA,F.h,CtrlVar.LevelSetMinIceThickness+1);
                        plot(xf/1000,yf/1000,'m','LineWidth',2);
                        
                        fMeshLSF.CurrentAxes.Title.String=[Temp(:)',{"White: LSF zero-line. Magneta:Ice cliff front"}];
                        
                    end
                    
                    % Par.RelativeVelArrowSize=10 ;
                    % QuiverColorGHG(MUA.coordinates(:,1)/1000,MUA.coordinates(:,2)/1000,F.ub,F.vb,Par) ;
                    
                    
                    %                         if contains(FileNameSubstring,'-1dIceShelf-')
                    %                             xlim([min(xc)-50e3  max(xc)+50e3]/1000)
                    %                             ylim([min(y) max(y)]/1e3) ;
                    %                         else
                    %                             if ~isnan(AxisLimits) ; axis(AxisLimits) ; end
                    %                         end
                    
                    xlabel('x (km)') ; ylabel('y (km)') ;
                    drawnow
                end
                
            case "-1dIceShelf-"
                %%
                
                
                
                [s,b,u,x]=AnalyticalOneDimentionalIceShelf(CtrlVar,MUA);
                yProfile=0 ;
                
                FigureName='flowline';
                FigFL=FindOrCreateFigure(FigureName) ;
                clf(FigFL)
                FigFL.InnerPosition=[100 700 939 665];
                hold on
                % point selection
                Iy=abs(MUA.coordinates(:,2)-yProfile)< 10000 ;
                
                xProfile=MUA.coordinates(Iy,1) ;
                [xProfile,Ix]=sort(xProfile) ;
                
                sProfile=F.s(Iy);
                bProfile=F.b(Iy);
                uProfile=F.ub(Iy) ;
                
                
                
                uProfile=uProfile(Ix) ;
                sProfile=sProfile(Ix);
                bProfile=bProfile(Ix);
                
                if isfield(F,'c') && ~isempty(F.c)
                    cProfile=F.c(Iy);
                    cProfile=cProfile(Ix);
                else
                    cProfile=[];
                end
                
                if isfield(F,'LSF') && ~isempty(F.LSF)
                    LSFProfile=F.LSF(Iy);
                    LSFProfile=LSFProfile(Ix);
                else
                    LSFProfile=[];
                end
                
                yyaxis left
                plot(xProfile/1000,sProfile,'bo-','DisplayName','$s$')
                hold on
                plot(xProfile/1000,bProfile,'go-','DisplayName','$b$')
                
                if contains(CtrlVar.Experiment,'analytical')
                    plot(x/1000,s,'b-','LineWidth',2,'DisplayName','$s$ analytical')
                    plot(x/1000,b,'g-','LineWidth',2,'DisplayName','$b$ analytical')
                end
                
                LSFmean=NaN;
                if ~isempty(LSFProfile)
                    
                    if isempty(LSFScale)
                        LSFScale=(max(sProfile)-min(bProfile))/(max(LSFProfile)-min(LSFProfile)) ;
                    end
                    
                    plot(xProfile/1000,LSFProfile*LSFScale,'m.-','DisplayName','$\varphi$ (scaled)')
                    
                    
                    if ~isempty(F.LSF)
                        [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F.LSF,0);
                        LSFmax=max(xc) ;
                        LSFmin=min(xc) ;
                        LSFmean=mean(xc) ;
                    end
                    
                    
                end
                
                ylabel('$z$ (m)','interpreter','latex')
                
                yyaxis right
                
                plot(xProfile/1000,uProfile,'ro-','DisplayName','$u$')
                hold on
                if contains(CtrlVar.Experiment,'analytical')
                    plot(x/1000,u,'r-','LineWidth',2,'DisplayName','$u$ analytical')
                end
                
                if ~isempty(cProfile)
                    plot(xProfile/1000,cProfile,'-k','DisplayName','$c$')
                end
                
                
                ylabel('$u$ (m/a)','interpreter','latex')
                
                title(sprintf('Profile along the medial line at t=%g',CtrlVar.time))
                xlabel('$x$ (km)','interpreter','latex') ;
                legend('interpreter','latex','Location','SouthEast')
                xlim([min(x)/1000 max(x)/1000]);
                
                ytt=ylim;
                
                plot([LSFmean,LSFmean]/1000,[ytt(1) ytt(2)],'k--')
                
                hold off
                
                
            case {"-mesh-speed-s-ab-","-mesh-speed-calving-level set-"}
                
                f4=FindOrCreateFigure("-mesh-speed-s-ab-",PlotScreenPosition);
                clf(f4)
                hold off
                
                
                
                % SP=tight_subplot(nPx,nPy) ;
                subplot(nPx,nPy,1)
                
                
                
                hold off
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.h);
                
                hold on
                PlotMuaMesh(CtrlVar,MUA,[],'w');
                
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r','LineWidth',2);
                if contains(PlotType,"-calving-")
                    [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,'w','LineWidth',2) ;
                end
                xlabel('x (km)') ; ylabel('y (km)') ;
                axis equal tight
                
                if ~isnan(AxisLimits) ; axis(AxisLimits) ; end
                title(sprintf('Ice thickness at t=%-g (yr)  #Ele=%-i, #Nodes=%-i, #nod=%-i',time,MUA.Nele,MUA.Nnodes,MUA.nod))
                title(cbar,'(m)')
                ax = gca;
                outerpos = ax.OuterPosition;
                ti = ax.TightInset;
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                ax.Position = [left bottom ax_width ax_height];
                
                hold off
                
                
                subplot(nPx,nPy,2);
                
                
                speed=sqrt(F.ub.*F.ub+F.vb.*F.vb);
                hold off
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,speed); title(sprintf('speed at t=%-g',time))
                %QuiverColorGHG(MUA.coordinates(:,1),MUA.coordinates(:,2),ub,vb,CtrlVar);
                hold on
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r','LineWidth',2);
                xlabel('x (km)') ; ylabel('y (km)') ; title(cbar,'(m/yr)')
                axis equal tight;
                if ~isnan(AxisLimits) ; axis(AxisLimits) ; end
                ax = gca;
                outerpos = ax.OuterPosition;
                ti = ax.TightInset;
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                ax.Position = [left bottom ax_width ax_height];
                if contains(PlotType,"-calving-")
                    [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,'w','LineWidth',2) ;
                end
                hold off
                
                subplot(nPx,nPy,3);
                
                
                hold off
                
                if contains(PlotType,"-calving-")
                    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.c);
                    title(sprintf('Calving Rate Field at t=%g',time))
                    hold on
                    [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,'w','LineWidth',2) ;
                else
                    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.s);
                    title(sprintf('surface at t=%-g',time))
                end
                
                hold on
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r','LineWidth',2);
                xlabel('x (km)') ; ylabel('y (km)') ;
                if ~isempty(cbar)
                    title(cbar,'(m/yr)')
                end
                axis equal tight ;
                if ~isnan(AxisLimits) ; axis(AxisLimits) ; end
                ax = gca;
                outerpos = ax.OuterPosition;
                ti = ax.TightInset;
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                ax.Position = [left bottom ax_width ax_height];
                
                subplot(nPx,nPy,4);
                
                
                hold off
                
                if contains(PlotType,"-level set-")
                    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.LSF/1000);
                    title(sprintf('Level Set Field at t=%-g',time))
                    hold on
                    [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,'w','LineWidth',2) ;
                else
                    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.ab);
                    title(sprintf('Basal melt at t=%-g',time))
                end
                
                hold on
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r','LineWidth',2);
                xlabel('x (km)') ; ylabel('y (km)') ;
                if ~isempty(cbar)
                    title(cbar,'(km)')
                end
                axis equal tight ;
                if ~isnan(AxisLimits) ; axis(AxisLimits) ; end
                hold off
                
                ax = gca;
                outerpos = ax.OuterPosition;
                ti = ax.TightInset;
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                ax.Position = [left bottom ax_width ax_height];
                
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
                caxis([-10 10])
                hold on ;
                [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'k');
                
                title(sprintf('dhdt at t=%-g (yr)',CtrlVar.time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
                title(colorbar,'(m)')
                if PlotMinThickLocations
                    plot(MUA.coordinates(ih,1)/CtrlVar.PlotXYscale,MUA.coordinates(ih,2)/CtrlVar.PlotXYscale,'.r');
                end
                
                
                
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
                
                
                
                
        end
        %
        %             if PlotType~="-1dIceShelf-"
        %                 axis equal tight ; axis(AxisLimits) ;
        %             end
        
        if CreateVideo
            iFrame=iFrame+1;
            if iFrame>1
                Frame = getframe(gcf);
                %Frame = hardcopy(hFig, '-opengl', '-r0');
                writeVideo(vidObj,Frame);
                hold off
            end
            
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

