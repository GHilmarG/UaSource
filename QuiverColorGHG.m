function [cbar,uvPlotScale,QuiverHandel]=QuiverColorGHG(x,y,u,v,CtrlVar,varargin)


% [cbar,uvPlotScale]=QuiverColorGHG(x,y,u,v,CtrlVar,varargin)
% a simple wrapper around quiver to generate coloured arrow field with a colorbar
%
% QuiverColorGHG(x,y,u,v)
% QuiverColorGHG(x,t,u,v,CtrlVar,varargin)
%
% CtrlVar.RelativeVelArrowSize                   : scaling factor for arrrow size, default value is 1
% CtrlVar.VelArrowColorSteps                     : number of coloring steps, default is 20
% CtrlVar.VelColorBarTitle                       : default value is '(m a^{-1})' ;
% CtrlVar.PlotXYscale                            : default value is 1
% CtrlVar.VelColorMap                            : default value is 'jet'
% CtrlVar.MinSpeedWhenPlottingVelArrows          : where speed is less, speed is set to this value, default value is zero
% CtrlVar.MinSpeedToPlot                         : where speed is less, speed is not plotted, default value is zero
% CtrlVar.VelPlotIntervalSpacing='lin'|'log10'   : lin or log10 vel scale (not sure if the colorbar is always correct for the log10 option...)
% CtrlVar.MaxPlottedSpeed=max(speed(:));
% CtrlVar.MinPlottedSpeed=min(speed(:));
% CtrlVar.SpeedTickLabels                        ; numerical array of values
% CtrlVar.QuiverColorPowRange                    ; when using log10 velocity bar, this is the range of magnitudes shown in colobar. 
%                                                  Default is  CtrlVar.QuiverColorPowRange=4, i.e. the smallest color is that of spee
%                                                  10^4 smaller than the largest speed
%                                                  Setting, for example, CtrlVar.QuiverColorPowRange=2 often gives better indication of changesin velocities.
% varargin is passed on to quiver
%
% Examples:
% QuiverColorGHG(x,y,ub,vb);
%
% velocities on top of FE mesh:
% PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
% hold on
% QuiverColorGHG(x,y,ub,vb,CtrlVar);
%
% Plot all velocities with same length arrows, colour coding shows actual speed
% speed=sqrt(ub.*ub+vb.*vb);
% CtrlVar.MinSpeedWhenPlottingVelArrows=max(speed);
% QuiverColorGHG(x,y,ub,vb,CtrlVar);
%  
%

if numel(x) ==0
    return
end

if size(u,1)> 1 && size(u,2)>1
    
   if size(u)==size(v)

       [X,Y]=meshgrid(x,y) ;
       x=X(:) ; y=Y(:) ; u=u(:) ; v=v(:);
       clear X Y
   end
   
end

x=x(:) ; y=y(:) ; u=u(:) ; v=v(:);

speed=sqrt(u.*u+v.*v); % speed is never scaled, so I can use speed to color velocity field based on values


if nargin>4
    
    if ~isfield(CtrlVar,'uvPlotScale')
        CtrlVar.uvPlotScale=[];
    end
    if ~isfield(CtrlVar,'QuiverColorPowRange')
        CtrlVar.QuiverColorPowRange=3; 
    end
    
    if ~isfield(CtrlVar,'RelativeVelArrowSize')
        CtrlVar.RelativeVelArrowSize=1; % larger value makes vel arrows larger
    end
    
    if ~isfield(CtrlVar,'VelColorMap')
        CtrlVar.VelColorMap='jet';
    end
    
    if ~isfield(CtrlVar,'PlotXYscale')
        CtrlVar.PlotXYscale=1;
    end
    
    if ~isfield(CtrlVar,'VelArrowColorSteps')
        N=50;
    else
        N=CtrlVar.VelArrowColorSteps;
    end
    
    if ~isfield(CtrlVar,'VelColorBarTitle')
        CtrlVar.VelColorBarTitle='(m a^{-1})' ;
    end
    
    if ~isfield(CtrlVar,'MinSpeedWhenPlottingVelArrows')
        CtrlVar.MinSpeedWhenPlottingVelArrows=0;
    end
    
    if ~isfield(CtrlVar,'MaxPlottedSpeed')
        
        if all(speed==0)
            ticks=logticks(speed,CtrlVar.QuiverColorPowRange);
            CtrlVar.MaxPlottedSpeed=max(ticks);
        else
            CtrlVar.MaxPlottedSpeed=max(speed(:))*1.001;
        end
    end
    
    
    if ~isfield(CtrlVar,'VelPlotIntervalSpacing')
        CtrlVar.VelPlotIntervalSpacing='lin';
    end
    
    if ~isfield(CtrlVar,'MinPlottedSpeed')
        switch CtrlVar.VelPlotIntervalSpacing
            case 'log10'
                
                ticks=logticks(speed,CtrlVar.QuiverColorPowRange);
                CtrlVar.MinPlottedSpeed=min(ticks);
%                 temp1=floor(10*min(speed(:)))/10;
%                 temp2=CtrlVar.MaxPlottedSpeed/1e3;
%                 CtrlVar.MinPlottedSpeed=max([temp1 temp2]);
                
            case 'lin'
                CtrlVar.MinPlottedSpeed=min(speed(:));
        end
    end
    
    if ~isfield(CtrlVar,'MinSpeedToPlot')
        CtrlVar.MinSpeedToPlot=0;
    end
    
    if ~isfield(CtrlVar,'SpeedTickLabels')
        CtrlVar.SpeedTickLabels=[];
    end
else
    
    CtrlVar.RelativeVelArrowSize=1; % larger value makes vel arrows larger
    N=50;
    CtrlVar.uvPlotScale=[];
    CtrlVar.VelColorBarTitle='(m a^{-1})' ;
    CtrlVar.PlotXYscale=1;
    CtrlVar.VelColorMap='jet';
    CtrlVar.MinSpeedWhenPlottingVelArrows=0;
    CtrlVar.MaxPlottedSpeed=max(speed(:));
    CtrlVar.MinPlottedSpeed=min(speed(:));
    CtrlVar.VelPlotIntervalSpacing='lin';
    CtrlVar.MinSpeedToPlot=0;
    CtrlVar.SpeedTickLabels=[];
    CtrlVar.QuiverColorPowRange=3;
    
end


cmap=colormap(sprintf('%s(%i)',CtrlVar.VelColorMap,N));


switch CtrlVar.VelPlotIntervalSpacing
    
    case 'log10'
        sp=logspace(log10(CtrlVar.MinPlottedSpeed),log10(CtrlVar.MaxPlottedSpeed),N+1);
    case 'lin'
        sp=linspace(CtrlVar.MinPlottedSpeed,CtrlVar.MaxPlottedSpeed,N+1); sp(1)=sp(1);
    otherwise
        fprintf(' which case {log10,lin}?' )
        error('QuiverColorGHG:VelPlotIntervalSpacing','case not reckognized')
end

% scaling of velocity to get resonably sized arrows

Np=(sqrt(numel(x))/10);
if Np<30 ; Np=30 ; end

ps=min([max(x)-min(x) max(y)-min(y)])/CtrlVar.PlotXYscale*CtrlVar.RelativeVelArrowSize/Np;

% the maximum velocity, when plotted, is multiplied by ps
%  uvPlotScale=max(speed)/ps  ;

%    if strcmp(CtrlVar.VelPlotIntervalSpacing,'log10')==1
%
%
%        uvPlotScale= speed/ps  ;
%
%
%    else
if isempty(CtrlVar.uvPlotScale)
    uvPlotScale= CtrlVar.MaxPlottedSpeed/ps  ;
else
    uvPlotScale=CtrlVar.uvPlotScale;
end
u=u/uvPlotScale; v=v/uvPlotScale;
CtrlVar.MinSpeedWhenPlottingVelArrows=CtrlVar.MinSpeedWhenPlottingVelArrows/uvPlotScale;

% end

for J=1:N
    if J==1
        I=speed <= sp(J+1) & speed>CtrlVar.MinSpeedToPlot;
    else
        I=speed>sp(J) & speed <= sp(J+1) & speed>CtrlVar.MinSpeedToPlot;
    end
    
    uplot=u(I) ; vplot=v(I);
    % where speed is less than CtrlVar.MinPlottingSpeed;
    % speed is set equal to CtrlVar.MinPlottingSpeed;
    speedplot=sqrt(u(I).*u(I)+v(I).*v(I));
    Ind=speedplot < CtrlVar.MinSpeedWhenPlottingVelArrows;
    uplot(Ind)=uplot(Ind)*CtrlVar.MinSpeedWhenPlottingVelArrows./speedplot(Ind);
    vplot(Ind)=vplot(Ind)*CtrlVar.MinSpeedWhenPlottingVelArrows./speedplot(Ind);
    
    QuiverHandel=quiver(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,uplot,vplot,0,'color',cmap(J,:),varargin{:}) ; hold on
end


nPowRange=CtrlVar.QuiverColorPowRange;

if verLessThan('matlab','8.4')
    
    %pre 2014b version
    cbar=colorbar; title(cbar,CtrlVar.VelColorBarTitle)   ;
    
    if strcmp(CtrlVar.VelPlotIntervalSpacing,'log10')==1
        
        ms=fix(log10(CtrlVar.MaxPlottedSpeed));
        ticklabel=logspace(0,ms,ms+1);
        tickpos=1+N*log10(ticklabel)/log10(sp(N+1));
        set(cbar,'Ytick',tickpos,'YTicklabel',ticklabel);
    else
        caxis([0 max(sp)]);
    end
    
    
else
    
    if strcmp(CtrlVar.VelPlotIntervalSpacing,'log10')==1
        
        
        if ~isempty(CtrlVar.SpeedTickLabels)
            ticklabel=log10(CtrlVar.SpeedTickLabels);
        else
            % this is for more pleasing interval between labels (but needs to be improved)
            %D=(max(sp)-min(sp))/10 ; D=10.^round(log10(D)) ; ticklabel=unique(D*round(sp/D));
            ticklabel=logticks(speed,nPowRange);
        end
        
        
        tickpos=(log10(ticklabel)-min(log10(sp)))/(max(log10(sp))-min(log10(sp)));
        
        
    else
        
        if ~isempty(CtrlVar.SpeedTickLabels)

            ticklabel=CtrlVar.SpeedTickLabels;
        else
            
            % this is for more pleasing interval between labels
            D=(max(sp)-min(sp));
            if D==0 ;  % special case if all values same or all values zero
                D=mean(sp);
                if D==0
                    ticklabel=[-1 0 1];
                    tickpos=(ticklabel+1)/2;
                else
                    ticklabel=[D/2 D 1.5*D];
                    tickpos=(ticklabel-D/2)/D;
                end
            else
                
                D=double(10.^floor(log10(D)))/2 ;
                
                first=D*floor(sp(1)/D);
                last=D*ceil(sp(end)/D);
                
                ticklabel=first:D:last;
                tickpos=(ticklabel-min(sp))/(max(sp)-min(sp));
            end
            %ticklabel=unique(D*round(sp/D));
            
        end
        
        
    end
    
    if any(isnan(tickpos))
        warning('QuiverColorGHG:NANinTickPos','calculated positions of ticks on the colorbar contain NaNs')
    end
    
    [tickpos,ia]=unique(tickpos); 
    ticklabel=ticklabel(ia);
    
    cbar=colorbar ; 
    cbar.TickLabels=ticklabel ; 

    if ~any(isnan(tickpos))
        cbar.Ticks=tickpos;
        title(cbar,CtrlVar.VelColorBarTitle)   ;
    end
    
    axis equal
    
    axis([min(x)/CtrlVar.PlotXYscale max(x)/CtrlVar.PlotXYscale min(y)/CtrlVar.PlotXYscale max(y)/CtrlVar.PlotXYscale])

end

end

