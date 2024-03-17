function [cbar,QuiverHandel,Par]=QuiverColorGHG(x,y,u,v,Par,varargin)

%% Plot velocity using colours/colors
%
% Just a wrapper around quiver to generate coloured arrow field with a colorbar.
%
%   [cbar,QuiverHandel,Par]=QuiverColorGHG(x,y,u,v,Par,varargin)
%
% x , y , u , v : vectors of same length. But if using regular velocity grid, x and y can
% be grid vectors. (In matlab speak grid vectors xg and yg are a set of vectors that serve
% as a compact representation of a grid in ndgrid format. For example, [X,Y] =
% ndgrid(xg,yg) )
%
%
%   Par.RelativeVelArrowSize                   : affects the size of the velocity arrows. 
%                                                by default Par.RelativeVelArrowSize=1.
%                                                Increase value for larger arrows.
%   Par.QuiverColorSpeedLimits=[min max]        : Speed range being colored, leave empty for auto.
%                                                Note however that when plotting
%                                                using log10 scaling the min speed
%                                                colored will be some fraction of the
%                                                max speed based on the value of
%                                              Par.QuiverColorPowRange
%                                            : scaling factor for arrow size, default value is 1
% Par.VelArrowColorSteps                     : number of coloring steps, default is 20
% Par.VelColorBarTitle                       : default value is '(m a^{-1})' ;
% Par.PlotXYscale                            : default value is 1
% Par.VelColorMap                            : default value is 'jet'
% Par.MinSpeedToPlot                         : where speed is less, speed is not plotted, default value is zero
% Par.VelPlotIntervalSpacing='lin'|'log10'   : lin or log10 vel scale
% Par.MaxPlottedSpeed                        : When plotting, speed above this value is set equal to this value, i.e. this is the maximum plotted speed
%                                                  Default is max(speed(:))
% Par.MinPlottedSpeed                        : When plotting, speed below this value is set equal to this value, i.e. this is the mainimum plotted speed
%                                                  Default is min(speed(:)).
%                                                  However, if using log10 the minimum plotted speed is never smaller than 10^QshouldiverColorPowRange times MaxPlottedSpeed
% Par.SpeedTickLabels                        : numerical array of values
% Par.QuiverColorPowRange                    : when using log10 velocity bar, this is the greates possible range of magnitudes shown in colobar.
%                                              Default is
%                                              Par.QuiverColorPowRange=3, i.e.
%                                              the smallest colored speed is 
%                                              10^3 smaller than the largest speed
% Par.QuiverSameVelocityScalingsAsBefore      : set to true (ie 1) to get same velocity scalings as in previous call. 
%                                           
%
% varargin is passed on to quiver
%
% *Examples:*
%
%  figure
%  load('CrackRestartfileExample.mat','CtrlVarInRestartFile','MUA','F','BCs','GF')
%  QuiverColorGHG(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub,F.vb);
%
% Plot all velocites with arrows of equal length:
%
%  load('CrackRestartfileExample.mat','CtrlVarInRestartFile','MUA','F','BCs','GF')
%   CtrlVar=CtrlVarInRestartFile;
%  speed=sqrt(F.ub.*F.ub+F.vb.*F.vb);
%  Par.MinPlottedSpeed=max(speed);
%  Par.VelColorBarTitle=' ';
%  Par.PlotXYscale=CtrlVar.PlotXYscale
%  figure
%  QuiverColorGHG(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub,F.vb,Par);
%  hold on ; PlotMuaBoundary(CtrlVar,MUA,'k')
%
% Plot velocities at approx equally spaced intervals:
%
%   load('CrackRestartfileExample.mat','CtrlVarInRestartFile','MUA','F','BCs','GF')
%   CtrlVar=CtrlVarInRestartFile;
%   x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2); 
%   [X,Y]=ndgrid(linspace(min(x),max(x),20),linspace(min(y),max(y),20));
%   I=nearestNeighbor(MUA.TR,[X(:) Y(:)]);  % finds nodes within computational grid closest to the regularly scape X and Y grid points.
%   FigVelocities=figure; 
%   Par.PlotXYscale=CtrlVar.PlotXYscale ; 
%   Par.MinPlottedSpeed=0; 
%   QuiverColorGHG(MUA.coordinates(I,1),MUA.coordinates(I,2),F.ub(I),F.vb(I),Par);
%   hold on ; PlotMuaBoundary(CtrlVar,MUA,'k')
%   
%
% Plot velocities using logarithmic scaling.
%
%   load('CrackRestartfileExample.mat','CtrlVarInRestartFile','MUA','F','BCs','GF')
%   CtrlVar=CtrlVarInRestartFile;
%   FigVelocities=figure; 
%   CtrlVar.VelPlotIntervalSpacing='log10' ; 
%   QuiverColorGHG(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub,F.vb,CtrlVar);
%   hold on ; PlotMuaBoundary(CtrlVar,MUA,'k')
%
% Plotting velocities on top of FE mesh:
% 
%   load('CrackRestartfileExample.mat','CtrlVarInRestartFile','MUA','F','BCs','GF')
%   CtrlVar=CtrlVarInRestartFile;
%   PlotMuaMesh(CtrlVar,MUA)
%   hold on
%   QuiverColorGHG(x,y,F.ub,F.vb,CtrlVar);
%
%
%
% Two calls with same velocity scaling:
%
%   [cbar,~,Par]=QuiverColorGHG(x,y,u,v,Par)  ; % first call, here Par is not strickly needed as an input
%   Par.QuiverSameVelocityScalingsAsBefore=1;
%   QuiverColorGHG(x,y,u,v,Par) ; % second call uses same scalings as previous one
%                                        
%
% Note: When doing further contour plots on top of velocity plot, matlab will possibly change the
% limits of the colorbar and the position of the ticklables will no longer be correct.
% If this happens then reset range and ticks, for example:
%
%   [~,~,Par]=QuiverColorGHG(x,y,ub,vb);
% ...some other plots that affect the colorbar labels, colors, etc...
%   cbar=colorbar;
%   cbar.Ticks=Par.QuiverTicks*(cbar.Limits(2)-cbar.Limits(1))+cbar.Limits(1);
%   cbar.TickLabels=Par.QuiverTickLabels;
%   title(cbar,'(m/d)')   ;
%%

persistent SpeedPlotIntervals uvPlotScale QuiverTickLabels QuiverTicks QuiverCmap



cbar=[];
QuiverHandel=[];

if isempty(x) || isempty(y) || isempty(u)   || isempty(v)
    return
end


if size(u,2)== 1 && (numel(x) ~= numel(y))
    error('Ua:QuiverColorGHG:xyDimensionsNotCompatible','x and y must have the same number of elements.')
end


if numel(u) ~= numel(v)
    error('Ua:QuiverColorGHG:uvDimensionsNotCompatible','u and v must have the same number of elements.')
end

%
% The expected typical useage is to plot one-dimensional arrays of velocites
% But u, v, x ,and y can also be given on a grid.
% If u and v is given on a grid, create vectors



if size(u,1)> 1 && size(u,2)>1 && ((size(x,2)==1 && size(y,2)==1)  || (size(x,1)==1 && size(y,1)==1 ))
    
    % u and v are matrices , x and y are vectors
    x=x(:) ; y=y(:) ;
    
    if size(u)==size(v)
        
        if (size(x,1)==size(u,1)) && (size(y,1)==size(u,2))
            
            [X,Y]=ndgrid(x,y) ;
            x=X(:) ; y=Y(:) ; u=u(:) ; v=v(:);
            clear X Y
            
        elseif (size(x,1)==size(u,2)) && (size(y,1)==size(u,1))
            
            warning('QuiverColorGHG:wrongdimensions','x and y are not grid vectors')
            
            [X,Y]=meshgrid(x,y) ;
            x=X(:) ; y=Y(:) ; u=u(:) ; v=v(:);
            clear X Y
            
        end
    end
end

x=x(:) ; y=y(:) ; u=u(:) ; v=v(:);

speed=sqrt(u.*u+v.*v); % speed is never scaled, so I can use speed to color velocity field based on values

if isinf(max(speed))
   
    fprintf(' Max of speed is infinite\n')
    
    I=isinf(speed);
    u(I)=NaN ; v(I)=NaN ; speed(I)=NaN;
    fprintf(' Seeting all such velocity values to NaN \n')
    
end


% now check Par fields, use all user-defined values where available, use default values otherwise
% and put in some reasonable for the remaining fields

if nargin>4 && ~isempty(Par)
    
    if ~isfield(Par,'QuiverColorPowRange')
        Par.QuiverColorPowRange=3;
    end
    
    if ~isfield(Par,'MaxPlottedSpeed')  || isempty(Par.MaxPlottedSpeed)
        
        if isfield(Par,'QuiverColorSpeedLimits')  && ~isempty(Par.QuiverColorSpeedLimits)
            Par.MaxPlottedSpeed=Par.QuiverColorSpeedLimits(2);
        else
            
            
            if all(speed==0)
                ticks=logticks(speed,Par.QuiverColorPowRange);
                Par.MaxPlottedSpeed=max(ticks);
                
            else
                %Par.MaxPlottedSpeed=max(speed(:))*1.001;
                Par.MaxPlottedSpeed=max(speed(:));
            end
        end
    end
    
    
    if ~isfield(Par,'MinPlottedSpeed')  || isempty(Par.MinPlottedSpeed)
        Par.MinPlottedSpeed=0;
    end
    
    if ~isfield(Par,'VelPlotIntervalSpacing') || isempty(Par.VelPlotIntervalSpacing)
        Par.VelPlotIntervalSpacing='lin';
    end
    
    
    if ~isfield(Par,'QuiverColorSpeedLimits')  || isempty(Par.QuiverColorSpeedLimits)
        
        switch Par.VelPlotIntervalSpacing
            case 'log10'
                
                ticks=logticks(speed,Par.QuiverColorPowRange,12);
                Par.QuiverColorSpeedLimits(1)=min(ticks);
                Par.QuiverColorSpeedLimits(2)=max(ticks);
            case 'lin'
                
                Par.QuiverColorSpeedLimits(1)=min(speed);
                Par.QuiverColorSpeedLimits(2)=max(speed);
        end
    end
    
    if ~isfield(Par,'uvPlotScale')
        Par.uvPlotScale=[];
    end
    
    
    if ~isfield(Par,'RelativeVelArrowSize') || isempty(Par.RelativeVelArrowSize)
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
        Par.VelColorBarTitle="($\mathrm{m \, yr^{-1}}$)" ; 
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
        Par.QuiverSameVelocityScalingsAsBefore=false ;
    end

    if Par.QuiverSameVelocityScalingsAsBefore

        if ~isempty(uvPlotScale)

            Par.uvPlotScale=uvPlotScale ;
            Par.SpeedPlotIntervals=SpeedPlotIntervals ;
            Par.QuiverTickLabels=QuiverTickLabels;
            Par.QuiverTicks=QuiverTicks ;
            Par.QuiverCmap=QuiverCmap;

        else
            Par.QuiverSameVelocityScalingsAsBefore=false ;
            fprintf('In QuiverColorGHG same scaling as in a previous call is requested, but the velocity scaling factor is not defined.\n')
            fprintf('Therefore, setting QuiverSameVelocityScalingsAsBefore=false \n')

        end
    end

else
    
    Par.RelativeVelArrowSize=1; % larger value makes vel arrows larger
    Par.VelArrowColorSteps=50;
    Par.uvPlotScale=[];
    Par.VelColorBarTitle="($\mathrm{m \, yr^{-1}}$)" ; 
    Par.PlotXYscale=1;
    Par.VelColorMap='jet';
    Par.MaxPlottedSpeed=max(speed(:));
    Par.MinPlottedSpeed=min(speed(:));
    Par.VelPlotIntervalSpacing='lin';
    Par.QuiverColorSpeedLimits(1)=min(speed(:));
    Par.QuiverColorSpeedLimits(2)=max(speed(:));
    Par.MinSpeedToPlot=0;
    Par.SpeedTickLabels=[];
    Par.QuiverColorPowRange=3;
    Par.QuiverCmap=[];
    Par.QuiverSameVelocityScalingsAsBefore=0;
    
end

%% 
% now all input variables should be OK
N=Par.VelArrowColorSteps;

if Par.QuiverSameVelocityScalingsAsBefore
    
    colormap(QuiverCmap)
    
    if isempty(Par.uvPlotScale) || ~isfield(Par,'uvPlotScale')
        fprintf('In QuiverColorGHG same scaling as in a previous call is requested, but the velocity scaling factor is not defined.\n')
        fprintf('Use parameters from previous call in this call!\n')
        error('Ua:QuiverColorGHG','The field uvPlotScale is not defiend')
    end
else
    





    if strcmp(Par.VelPlotIntervalSpacing,'log10')==1
        % create a `logarithmic' colormap

        if ~isnumeric(Par.VelColorMap)
            NN=10*N ;
            cmap=colormap(sprintf('%s(%i)',Par.VelColorMap,NN));
        else
            cmap=Par.VelColorMap ;
            NN=size(cmap,1);
        end
        index=fix((NN-1)*(exp((0:N-1)/(N-1))-1)/(exp(1)-1)+1);
        cmap=colormap(cmap(index,:));
    else

        if ~isnumeric(Par.VelColorMap)
            cmap=colormap(sprintf('%s(%i)',Par.VelColorMap,N));
        else
            cmap=Par.VelColorMap ;
            N=size(cmap,1);
        end
    end


    switch Par.VelPlotIntervalSpacing

        case 'log10'

            ticks=logticks(speed,Par.QuiverColorPowRange,12,Par.QuiverColorSpeedLimits);

            MinTick=min(ticks);

            if Par.QuiverColorSpeedLimits(1)<MinTick

                %Par.MinPlottedSpeed=MinTick;
                Par.QuiverColorSpeedLimits(1)=MinTick;
            end


            Par.SpeedPlotIntervals=logspace(log10(Par.QuiverColorSpeedLimits(1)),log10(Par.QuiverColorSpeedLimits(2)),N+1);

        case 'lin'

            %Par.SpeedPlotIntervals=linspace(Par.MinPlottedSpeed,Par.MaxPlottedSpeed,N+1);
            Par.SpeedPlotIntervals=linspace(Par.QuiverColorSpeedLimits(1),Par.QuiverColorSpeedLimits(2),N+1);

        otherwise
            fprintf(' which case {log10,lin}?' )
            error('QuiverColorGHG:VelPlotIntervalSpacing','case not reckognized')
    end




    Par.QuiverCmap=cmap;

    
    % scaling of velocity to get resonably sized arrows
    
    
    Np=(sqrt(numel(x))/10);
    
    if Np<30 ; Np=30 ; end
    
    ps=min([max(x)-min(x) max(y)-min(y)])/Par.PlotXYscale*Par.RelativeVelArrowSize/Np;
    Par.uvPlotScale = Par.QuiverColorSpeedLimits(2)/ps  ;
    
    
end

uplot=u ; vplot=v;

% set velocities within plotting range
Ind=speed > Par.MaxPlottedSpeed;
uplot(Ind)=uplot(Ind)*Par.MaxPlottedSpeed./speed(Ind);
vplot(Ind)=vplot(Ind)*Par.MaxPlottedSpeed./speed(Ind);

Ind=speed < Par.MinPlottedSpeed;
uplot(Ind)=uplot(Ind)*Par.MinPlottedSpeed./speed(Ind);
vplot(Ind)=vplot(Ind)*Par.MinPlottedSpeed./speed(Ind);



%fprintf('QuiverColorGHG: uvPlotScale=%f \n',Par.uvPlotScale)

uplot=uplot/Par.uvPlotScale; vplot=vplot/Par.uvPlotScale;

%SpeedPlot=sqrt(uplot.*uplot+vplot.*vplot);
%sc=log(1+SpeedPlot)./x  ; I=isnan(sc) ; sc(I)=1;

% end

for J=1:numel(Par.SpeedPlotIntervals)-1
    
    switch J
        case 1
            I=speed <= Par.SpeedPlotIntervals(J+1) & speed>Par.MinSpeedToPlot;
        case N
            I=speed>Par.SpeedPlotIntervals(J) & speed>Par.MinSpeedToPlot;
        otherwise
            I=speed>Par.SpeedPlotIntervals(J) & speed <= Par.SpeedPlotIntervals(J+1) & speed>Par.MinSpeedToPlot;
    end
    
    if numel(x(I))>0  % This should not really be needed, but veving resulting figures in matlab fig format
                      % results in errors in Matlab2016a
        QuiverHandel=quiver(x(I)/Par.PlotXYscale,y(I)/Par.PlotXYscale,uplot(I),vplot(I),0,...
            'color',Par.QuiverCmap(J,:),varargin{:}) ; hold on
    end
end


nPowRange=Par.QuiverColorPowRange;

if ~Par.QuiverSameVelocityScalingsAsBefore

    if verLessThan('matlab','8.4')

        %pre 2014b version
        cbar=colorbar; title(cbar,Par.VelColorBarTitle,"interpreter","latex")   ;

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
                %ticks=logticks(speed,12,Par.QuiverColorPowRange);
            end


            tickpos=(log10(ticklabel)-min(log10(Par.SpeedPlotIntervals)))/(max(log10(Par.SpeedPlotIntervals))-min(log10(Par.SpeedPlotIntervals)));


        else

            if ~isempty(Par.SpeedTickLabels)

                ticklabel=Par.SpeedTickLabels;
            else

                % this is for (hopefully) a more pleasing interval between labels
                D=(max(Par.SpeedPlotIntervals)-min(Par.SpeedPlotIntervals));
                if D==0  % special case if all values same or all values zero
                    D=mean(Par.SpeedPlotIntervals);
                    if D==0
                        ticklabel=[-1 0 1];
                        tickpos=(ticklabel+1)/2;
                    else
                        ticklabel=[D/2 D 1.5*D];
                        tickpos=(ticklabel-D/2)/D;
                    end
                else

                    D=double(10.^floor(log10(D)))/4 ;

                    first=D*floor(Par.SpeedPlotIntervals(1)/D);
                    last=D*ceil(Par.SpeedPlotIntervals(end)/D);

                    ticklabel=first:D:last;
                    tickpos=(ticklabel-min(Par.SpeedPlotIntervals))/(max(Par.SpeedPlotIntervals)-min(Par.SpeedPlotIntervals));
                    %tickpos=(ticklabel-first)/(last-first);

                    if numel(ticklabel)>10

                        ticklabel=ticklabel(1:2:end) ;
                        tickpos=tickpos(1:2:end) ;

                    end

                end
                %ticklabel=unique(D*round(sp/D));

            end


        end


        if any(isnan(tickpos))
            warning('QuiverColorGHG:NANinTickPos','calculated positions of ticks on the colorbar contain NaNs')
        end

        [tickpos,ia]=unique(tickpos);

        ticklabel=ticklabel(ia);
        colormap(Par.QuiverCmap)
        cbar=colorbar ;
        
        %cbar.TickLabels=ticklabel ;


        if ~any(isnan(tickpos))
            Ticks=tickpos*(cbar.Limits(2)-cbar.Limits(1))+cbar.Limits(1);
            %cbar.Ticks=Ticks; % tickpos*(cbar.Limits(2)-cbar.Limits(1))+cbar.Limits(1);
        end

        axis equal



        Par.QuiverTickLabels=ticklabel;
        Par.QuiverTicks=Ticks;

    end
end

if isempty(cbar)
    cbar=colorbar;
end
title(cbar,Par.VelColorBarTitle,"interpreter","latex")   ;
cbar.TickLabels=Par.QuiverTickLabels;
cbar.Ticks=Par.QuiverTicks;

% Par.QuiverColorSpeedLimits=[min(Par.SpeedPlotIntervals) max(Par.SpeedPlotIntervals)];
% Par.QuiverColorSpeedLimits=[];  % don't reuse these setting in next call


axis equal
% axis([min(x)/Par.PlotXYscale max(x)/Par.PlotXYscale min(y)/Par.PlotXYscale max(y)/Par.PlotXYscale])

uvPlotScale=Par.uvPlotScale ;
SpeedPlotIntervals=Par.SpeedPlotIntervals ;
QuiverTickLabels=Par.QuiverTickLabels;
QuiverTicks=Par.QuiverTicks ;
QuiverCmap=Par.QuiverCmap;


end

