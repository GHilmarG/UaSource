function [cbar,QuiverHandel,Par,Colorbar]=QuiverColorGHG(x,y,u,v,Par,varargin)

%% Plot velocity using colours/colors
%
% [cbar,uvPlotScale]=QuiverColorGHG(x,y,u,v,Par,varargin) a simple wrapper
% around quiver to generate coloured arrow field with a colorbar
%
% Example: QuiverColorGHG(x,y,u,v) QuiverColorGHG(x,t,u,v,Par,varargin)
%
% x , y , u , v : vectors of same length. But if using regular grid x and y can
% be grid vectors. (In matlab speak grid vectors are a set of vectors that serve
% as a compact representation of a grid in ndgrid format. For example, [X,Y] =
% ndgrid(xg,yg) returns a full grid in the matrices X and Y. You can represent
% the same grid using the grid vectors, xg and yg.)
%
%
% Par.RelativeVelArrowSize                   : scaling factor for arrrow size, default value is 1
% Par.VelArrowColorSteps                     : number of coloring steps, default is 20
% Par.VelColorBarTitle                       : default value is '(m a^{-1})' ;
% Par.PlotXYscale                            : default value is 1
% Par.VelColorMap                            : default value is 'jet'
% Par.MinSpeedToPlot                         : where speed is less, speed is not plotted, default value is zero
% Par.VelPlotIntervalSpacing='lin'|'log10'   : lin or log10 vel scale (not sure if the colorbar is always correct for the log10 option...)
%                                                :
% Par.MaxPlottedSpeed                        : When plotting speed above this value is set equal to this value, i.e. this is the maximum plotted speed
%                                                  Default is max(speed(:))
% Par.MinPlottedSpeed                        : When plotting speed below this value is set equal to this value, i.e. this is the mainimum plotted speed
%                                                  Default is min(speed(:)).
%                                                  However, if using log10 the minimum plotted speed is never smaller than 10^QuiverColorPowRange times MaxPlottedSpeed
% Par.SpeedTickLabels                        : numerical array of values
% Par.QuiverColorPowRange                    : when using log10 velocity bar, this is the greates possible range of magnitudes shown in colobar.
%                                                  Default is  Par.QuiverColorPowRange=3, i.e. the smallest color is that of spee
%                                                  10^3 smaller than the largest speed
%                                                  Setting, for example, Par.QuiverColorPowRange=2 narrows the plotted range
%
% Par.QuiverSameVelocityScalingsAsBefore      : set to true (ie 1) to get same velocity scalings as in previous call. In this case
%                                               Par from previous call must be given as an input.
%
% varargin is passed on to quiver
%
% Examples:
% QuiverColorGHG(x,y,ub,vb);
%
% velocities on top of FE mesh:
% PlotFEmesh(MUA.coordinates,MUA.connectivity,Par)
% hold on
% QuiverColorGHG(x,y,ub,vb,Par);
%
% Plot all velocities with same length arrows, colour coding shows actual speed
% speed=sqrt(ub.*ub+vb.*vb);
% Par.MinPlottedSpeed=max(speed);
% QuiverColorGHG(x,y,ub,vb,Par);
%
% Two calls with same velocity scaling:
% [c~,~,Par,Colorbar]=QuiverColorGHG(x,y,u,v,Par)  ; % first call, here Par is not needed as input
%  Par.QuiverSameVelocityScalingsAsBefore=1; 
% [c~,~,Par,Colorbar]=QuiverColorGHG(x,y,u,v,Par) ; % second call using same scalings as the previous one
%
% Note: When doing further contour plots on top of velocity plot, matlab will possibly change the
% limits of the colorbar and the position of the ticklables will no longer be correct.
% If this happens then reset range and ticks:
% cbar=colorbar;
% cbar.Ticks=Par.QuiverTicks*(cbar.Limits(2)-cbar.Limits(1))+cbar.Limits(1);
% cbar.TickLabels=Par.QuiverTickLabels;
% title(cbar,'(m/d)')   ;
%%


if numel(x) ==0
    return
end


%
% The expected typical useage is to plot one-dimentional arrays of velocites
% But u, v, x ,and y can also be given on a grid.
% If u and v is given on a grid, create vectors

if size(u,1)> 1 && size(u,2)>1
    
    if size(u)==size(v)
        
        %[X,Y]=meshgrid(x,y) ;
        [X,Y]=ndgrid(x,y) ;
        x=X(:) ; y=Y(:) ; u=u(:) ; v=v(:);
        clear X Y
    end
    
end

x=x(:) ; y=y(:) ; u=u(:) ; v=v(:);

speed=sqrt(u.*u+v.*v); % speed is never scaled, so I can use speed to color velocity field based on values


% now check Par fields, use all user-defined values where available, use default values otherwise
% and put in some reasonable for the remaining fields

if nargin>4 && ~isempty(Par)
    
    if ~isfield(Par,'QuiverColorPowRange')
        Par.QuiverColorPowRange=3;
    end
    
    if ~isfield(Par,'MaxPlottedSpeed')  || isempty(Par.MaxPlottedSpeed)
        
        if all(speed==0)
            ticks=logticks(speed,Par.QuiverColorPowRange);
            Par.MaxPlottedSpeed=max(ticks);
        else
            %Par.MaxPlottedSpeed=max(speed(:))*1.001;
            Par.MaxPlottedSpeed=max(speed(:));
        end
        
    end
    
    
    if ~isfield(Par,'VelPlotIntervalSpacing') || isempty(Par.VelPlotIntervalSpacing)
        Par.VelPlotIntervalSpacing='lin';
    end
    
    
    if ~isfield(Par,'MinPlottedSpeed')  || isempty(Par.MinPlottedSpeed)
        
        switch Par.VelPlotIntervalSpacing
            case 'log10'
                
                ticks=logticks(speed,12,Par.QuiverColorPowRange);
                Par.MinPlottedSpeed=min(ticks);
                %                 temp1=floor(10*min(speed(:)))/10;
                %                 temp2=Par.MaxPlottedSpeed/1e3;
                %                 Par.MinPlottedSpeed=max([temp1 temp2]);
                
            case 'lin'
                Par.MinPlottedSpeed=min(speed(:));
        end
    end
    
    if ~isfield(Par,'uvPlotScale')
        Par.uvPlotScale=[];
    end
    
    
    if ~isfield(Par,'RelativeVelArrowSize')
        Par.RelativeVelArrowSize=1; % larger value makes vel arrows larger
    end
    
    if ~isfield(Par,'VelColorMap')
        Par.VelColorMap='jet';
    end
    
    if ~isfield(Par,'PlotXYscale')
        Par.PlotXYscale=1;
    end
    
    if ~isfield(Par,'VelArrowColorSteps')
        Par.VelArrowColorSteps=50;
    end
    
    if ~isfield(Par,'VelColorBarTitle')
        Par.VelColorBarTitle='(m a^{-1})' ;
    end
    
    
    
    if ~isfield(Par,'MinSpeedToPlot')
        Par.MinSpeedToPlot=0;
    end
    
    if ~isfield(Par,'SpeedTickLabels')
        Par.SpeedTickLabels=[];
    end
    
    if ~isfield(Par,'QuiverCmap')
        Par.QuiverCmap=[];
    end
    
    if ~isfield(Par,'QuiverSameVelocityScalingsAsBefore')
        Par.QuiverSameVelocityScalingsAsBefore=0;
    end
    
else
    
    Par.RelativeVelArrowSize=1; % larger value makes vel arrows larger
    Par.VelArrowColorSteps=50;
    Par.uvPlotScale=[];
    Par.VelColorBarTitle='(m a^{-1})' ;
    Par.PlotXYscale=1;
    Par.VelColorMap='jet';
    Par.MaxPlottedSpeed=max(speed(:));
    Par.MinPlottedSpeed=min(speed(:));
    Par.VelPlotIntervalSpacing='lin';
    Par.MinSpeedToPlot=0;
    Par.SpeedTickLabels=[];
    Par.QuiverColorPowRange=3;
    Par.QuiverCmap=[];
    Par.QuiverSameVelocityScalingsAsBefore=0;
    
end

%% now all input variables should be OK
% 
N=Par.VelArrowColorSteps;

if Par.QuiverSameVelocityScalingsAsBefore
    
    colormap(Par.QuiverCmap)
    
else
    
    if ~isfield(Par,'SpeedPlotIntervals') || isempty(Par.SpeedPlotIntervals)
        
        switch Par.VelPlotIntervalSpacing
            
            case 'log10'
                
                ticks=logticks(speed,Par.QuiverColorPowRange);
                
                MinTick=min(ticks);
                
                if Par.MinPlottedSpeed<MinTick
                    
                    Par.MinPlottedSpeed=MinTick;
                    
                end
                
                
                Par.SpeedPlotIntervals=logspace(log10(Par.MinPlottedSpeed),log10(Par.MaxPlottedSpeed),N+1);
                
            case 'lin'
                Par.SpeedPlotIntervals=linspace(Par.MinPlottedSpeed,Par.MaxPlottedSpeed,N+1);
                
            otherwise
                fprintf(' which case {log10,lin}?' )
                error('QuiverColorGHG:VelPlotIntervalSpacing','case not reckognized')
        end
        
        
    end
    
    %%
    % Now all Par fields have been checked or set to some reasonable values
    %
    
    
    
    if strcmp(Par.VelPlotIntervalSpacing,'log10')==1
        % create a `logarithmic' colormap
        NN=10*N ;
        cmap=colormap(sprintf('%s(%i)',Par.VelColorMap,NN));
        index=fix((NN-1)*(exp([0:N-1]/(N-1))-1)/(exp(1)-1)+1);
        cmap=colormap(cmap(index,:));
    else
        cmap=colormap(sprintf('%s(%i)',Par.VelColorMap,N));
    end
    
    Par.QuiverCmap=cmap;
    
    
    % set velocities within plotting range
    Ind=speed > Par.MaxPlottedSpeed;
    u(Ind)=u(Ind)*Par.MaxPlottedSpeed./speed(Ind);
    v(Ind)=v(Ind)*Par.MaxPlottedSpeed./speed(Ind);
    
    Ind=speed < Par.MinPlottedSpeed;
    u(Ind)=u(Ind)*Par.MinPlottedSpeed./speed(Ind);
    v(Ind)=v(Ind)*Par.MinPlottedSpeed./speed(Ind);
    
    % scaling of velocity to get resonably sized arrows
    
    
    Np=(sqrt(numel(x))/10);
    
    if Np<30 ; Np=30 ; end
    
    ps=min([max(x)-min(x) max(y)-min(y)])/Par.PlotXYscale*Par.RelativeVelArrowSize/Np;
    Par.uvPlotScale = Par.MaxPlottedSpeed/ps  ;
    
    
end



fprintf('QuiverColorGHG: uvPlotScale=%f \n',Par.uvPlotScale)

uplot=u/Par.uvPlotScale; vplot=v/Par.uvPlotScale;

% end

for J=1:N
    
    switch J
        case 1
            I=speed <= Par.SpeedPlotIntervals(J+1) & speed>Par.MinSpeedToPlot;
        case N
            I=speed>Par.SpeedPlotIntervals(J) & speed>Par.MinSpeedToPlot;
        otherwise
            I=speed>Par.SpeedPlotIntervals(J) & speed <= Par.SpeedPlotIntervals(J+1) & speed>Par.MinSpeedToPlot;
    end
    
    QuiverHandel=quiver(x(I)/Par.PlotXYscale,y(I)/Par.PlotXYscale,uplot(I),vplot(I),0,...
        'color',Par.QuiverCmap(J,:),varargin{:}) ; hold on
end


nPowRange=Par.QuiverColorPowRange;

if verLessThan('matlab','8.4')
    
    %pre 2014b version
    cbar=colorbar; title(cbar,Par.VelColorBarTitle)   ;
    
    if strcmp(Par.VelPlotIntervalSpacing,'log10')==1
        
        ms=fix(log10(Par.MaxPlottedSpeed));
        ticklabel=logspace(0,ms,ms+1);
        tickpos=1+N*log10(ticklabel)/log10(Par.SpeedPlotIntervals(N+1));
        set(cbar,'Ytick',tickpos,'YTicklabel',ticklabel);
    else
        caxis([0 max(Par.SpeedPlotIntervals)]);
    end
    
    
else
    
    if strcmp(Par.VelPlotIntervalSpacing,'log10')==1
        
        
        if ~isempty(Par.SpeedTickLabels)
            ticklabel=log10(Par.SpeedTickLabels);
        else
            % this is for more pleasing interval between labels (but needs to be improved)
            %D=(max(sp)-min(sp))/10 ; D=10.^round(log10(D)) ; ticklabel=unique(D*round(sp/D));
            ticklabel=logticks(speed,nPowRange,12,Par.SpeedPlotIntervals);
        end
        
        
        tickpos=(log10(ticklabel)-min(log10(Par.SpeedPlotIntervals)))/(max(log10(Par.SpeedPlotIntervals))-min(log10(Par.SpeedPlotIntervals)));
        
        
    else
        
        if ~isempty(Par.SpeedTickLabels)
            
            ticklabel=Par.SpeedTickLabels;
        else
            
            % this is for (hopefully) a more pleasing interval between labels
            D=(max(Par.SpeedPlotIntervals)-min(Par.SpeedPlotIntervals));
            if D==0 ;  % special case if all values same or all values zero
                D=mean(Par.SpeedPlotIntervals);
                if D==0
                    ticklabel=[-1 0 1];
                    tickpos=(ticklabel+1)/2;
                else
                    ticklabel=[D/2 D 1.5*D];
                    tickpos=(ticklabel-D/2)/D;
                end
            else
                
                D=double(10.^floor(log10(D)))/2 ;
                
                first=D*floor(Par.SpeedPlotIntervals(1)/D);
                last=D*ceil(Par.SpeedPlotIntervals(end)/D);
                
                ticklabel=first:D:last;
                tickpos=(ticklabel-min(Par.SpeedPlotIntervals))/(max(Par.SpeedPlotIntervals)-min(Par.SpeedPlotIntervals));
                %tickpos=(ticklabel-first)/(last-first);
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
    %cbar.TickLabels=ticklabel ;
    
    
    if ~any(isnan(tickpos))
        Ticks=tickpos*(cbar.Limits(2)-cbar.Limits(1))+cbar.Limits(1);
        %cbar.Ticks=Ticks; % tickpos*(cbar.Limits(2)-cbar.Limits(1))+cbar.Limits(1);
    end
    
    
    
    axis equal
    
    axis([min(x)/Par.PlotXYscale max(x)/Par.PlotXYscale min(y)/Par.PlotXYscale max(y)/Par.PlotXYscale])
    
    Par.QuiverTickLabels=ticklabel;
    Par.QuiverTicks=Ticks;
    Par.QuiverLimits=cbar.Limits;
end


cbar=colorbar ;
title(cbar,Par.VelColorBarTitle)   ;
cbar.TickLabels=Par.QuiverTickLabels;
cbar.Ticks=Par.QuiverTicks;

Colorbar.handle=cbar;
Colorbar.Ticklabels=Par.QuiverTickLabels;
Colorbar.Ticks=Par.QuiverTicks;



end

