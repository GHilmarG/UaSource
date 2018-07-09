

%%
Experiment='ice1rr';

list=dir(['*',Experiment,'*.mat']);
nFiles=length(list);


vidObj = VideoWriter(Experiment,'MPEG-4');
vidObj.FrameRate=2;
%vidObj.FrameRate=1;   % frames per sec
open(vidObj);


TRI=[] ; DT=[] ; iFile=1; iFrame=0 ; dt=10;
AspectRatio=2;
ViewAndLight(1)=170 ;  ViewAndLight(2)= 60;
ViewAndLight(3)=30 ;  ViewAndLight(4)=50;
%nFiles=2;

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
        
    

    
    
   
    B=b;
    
    close all
    Fig=figure;
    
    %Fig.Position=[962 42 958 962];
    
    
    
    CtrlVar.ThicknessCutOffForPlotting=15; 
    ViewAndLight(1)=ViewAndLight(1)+0 ; 
    ViewAndLight(2)=ViewAndLight(2)+0 ; 
    hold off
    
    [TRI,DT]=Plot_sbB(CtrlVar,MUA,s,b,B,TRI,DT,AspectRatio,ViewAndLight);
    xlabel('(km)' ) 
    ylabel('(km)' )
    
    
    iFrame=iFrame+1;
    Frame = getframe(gcf);
    %Frame = hardcopy(hFig, '-opengl', '-r0');
    writeVideo(vidObj,Frame);
    hold off
    
    end
    
    iFile=iFile+1;
   
end


close(vidObj);
