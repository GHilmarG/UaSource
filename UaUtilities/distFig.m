function distFig(varargin)
% ===== Syntax ============================================================
% distFig(...,'Screen',Value)
% distFig(...,'Position',Value)
%
% distFig(...,'Rows',Value)
% distFig(...,'Columns',Value)
%
% distFig(...,'Not',Value)
% distFig(...,'Only',Value)
% distFig(...,'Offset',Value)
%
% distFig(...,'Adjust',Value)
% distFig(...,'Extra',Value)
% distFig(...,'Menu',Value)
% distFig(...,'Transpose',Value)
% distFig(...,'Simulink',Value)
% distFig(...,'Tight',Value)
% distFig(...,'Scale',Value)
% distFig(...,'ScaleOrigin',Value)
% distFig(...,'FreezeSize',Value)
% distFig(...,'Order',Value)
%
% ===== Description =======================================================
% distFig(...,'Screen',Value) assigns where the figures will be
% distributed.
% Synonyms for 'Screen' are 'Scr' and 'S'. Value can be:
% 'C' / 'Center' (Default)
% 'W' / 'West' / 'L' / 'Left'
% 'E' / 'East' / 'R' / 'Right'
% 'N' / 'North' / 'T' / 'Top'
% 'S' / 'South' / 'B' / 'Bottom'
% 'Main' / 'Primary'
% 'Ext' / 'External' / 'Secondary'
%
% distFig(...,'Position',Value) assigns in which part of the screen the
% figures will be distributed.
% Synonyms for 'Position' are 'Pos' and 'P'. Value can be:
% 'C' / 'Center' (Default)
% 'W' / 'West' / 'L' / 'Left'
% 'E' / 'East' / 'R' / 'Right'
% 'N' / 'North' / 'T' / 'Top'
% 'S' / 'South' / 'B' / 'Bottom'
% 'NW'
% 'NE'
% 'SW'
% 'SE'
%
% distFig(...,'Rows',Value) assigns how many rows the figures will be
% distributed on.
% Synonyms for 'Rows' are 'R'. Value can be:
% 1...n
% 'Auto' / -1
% 'Auto' indicates that it automatically calculates the number of required
% rows.
%
% distFig(...,'Columns',Value) assigns how many rows the figures will be
% distributed on.
% Synonyms for 'Columns' are 'Cols' and 'C'. Value can be:
% 1...n
% 'Auto' / -1
% 'Auto' indicates that it automatically calculates the number of required
% columns.
%
% distFig(...,'Not',Value) excludes specified figures from the distribution
% list. Value must be an array with the excluded figure numbers. The
% default value is []. The input argument 'gcf' (get current figure) can be
% used as well. 
%
% distFig(...,'Only',Value) does only distribute specified figures. Value
% must be an array with the figure which will be distributed. The
% default value is []. The input argument 'gcf' (get current figure) can be
% used as well. 
%
% distFig(...,'Offset',Value) can be used to shift all figures on the
% distribution list. Synonyms for 'Offset' is 'O'. Value must be an integer
% larger or equal to zero. The default offset value is 0.
%
% distFig(...,'Adjust',Value) can be used to adjust all figure positions. 
% For instance, if all figures should be moved 20 pixels to the right,
% which is not possible with the other features, Value should be [20,0].
% If Value is below 1, the scalars indicate the relative adjustment to the
% screen. If the selected monitor is [1920,1080] pixels, a Values of 
% [0.25,-0.5] moves all figures with [480,540] pixels.
% Synonyms for 'Ajust' is 'A'. Value can be:
% [Integer,Integer] (In order to adjust by a specific ammount)
% [Scalar,Scalar]( < 1) (To adjust a relative ammount of the monitor size)
% The default value is [0,0].
% 
% distFig(...,'Extra',Value) is used to handle the extra figures, which
% might not fit into 'Rows' * 'Cols'. For instance, 5 figure are to be
% distributed on a 2x2 = 4 area. The fifth figure is then places on top of
% the first figure if 'Extra' is 'restart'. It can be ignored by setting
% 'Extra' to 'ignore', in which the fifth figure will not be distributed.
%
% distFig(...,'Menu',Value) can be used to modify the menu bar.
% Value can be:
% 'figure'
% 'none'
%
% distFig(...,'Transpose',Value) can be used to tranpose the order of
% distributing. The values must be logical and the default value is false.
%
% distFig(...,'Simulink',Value) can be used to handle simulink figures.
% Value can be:
% 'include' (Default)
% 'exclude'
% 'only'
% If value is true, simulink figures will be distributed along with the
% other figures. If value if 'exclude' simulink figures will be
% excluded from the distribution list. If value is 'only', only the
% simulink figures will be distributed and not the ordinary figures.
%
% distFig(...,'Tight',Value) can be used to stretch the figure to the
% border of the figure, so no unnecessary "whitespace" is present.
% Value must be logical and the default value if false.
%
% distFig(...,'Scale',Value) can be used to scale the figure matrix. Both
% horizontal and vertical scaling can be applied, but also similar scaling.
% Value is a single scalar for similar scaling, but is a 2x1 matrix for
% specific scaling in each direction. For instance, distFig('Scale',0.5)
% scales all figures by 50% relative to the upper left corner of the 
% distribution area (Default Scale-origin). distFig('Scale',[0.8,0.4])
% scales the figures by 80% and 40% in the horizontal and vertical 
% directions respectivly relative to the default scale-origin. 
% 
% distFig(...,'ScaleOrigin',Value) can be used in combination with the 
% 'Scale' option. Value is a 2x1 matrix with the origin. If these values
% are between 0 and 1 the relative monitor size is used, whereas larger
% values should be rounded indicating pixel values. For instance, 
% distFig('Scale',0.5,'ScaleOrigin',[0.5,1]) scales the figures relative 
% to the midpoint on the right hand side of the selected monitor. For a
% monitor with resolution 1920x1080, this would be equal to 
% distFig('Scale',0.5,'ScaleOrigin',[1920,960]).
% The default scale origin is in the upper left corner of the distribution
% area.
% Synonyms for 'ScaleOrigin' are 'Origin', 'Orig' and 'Scaleorig'.
% 
% distFig(...,'FreezeSize',Value) can be used to freeze the sizes of the
% figures before distribution. Currently, this function anchors the figures
% to the upper left corner of the original position, yet keeping the 
% original size - it does not calculate new positions, so the figures may
% overlap. 
% Synonym for 'FreezeSize' is 'Freeze'. 
% 
% distFig(...,'Order',Value) can be used to control the plotting order. For
% instance, distFig('Only',[5,6,7,8],'Order',[2,1,4,3]) will plot figure
% 5,6,7 and 8 in the order of 6,5,8 and 7. The ordering argument is thus
% related to the input of 'Only'. 
% 
% ===== Examples ==========================================================
% distFig();
% This will distribute all open figures on the primary screen in two rows.
%
% distFig('Screen','Left','Position','Right','Only',[1,2,4])
% This will only distribute figure 1, 2 and 4 on the right part of the left
% screen.
%
% distFig('Offset',2,'Not',[1,2])
% This will distribute all figure but figure 1 and 2 in the same pattern as
% distFig(), but figure 1 and 2 will not be distributed, but instead there will
% be blank spots where they would have been distributed.
% 
% ===== Note ==============================================================
% If you find any errors or bugs with this version please contact me on
% AndersSSimonsen@GMail.com or write a comment on file exchange - then I'll
% take a look at it. :)
% 
% This function will only work in MATLAB 2014b and onwards due to the
% massive overhaul in the figure handle system - I will not support
% previous versions.

% =========================================================================
% ===== Version ===========================================================
% =========================================================================

% Anders Schou Simonsen, AndersSSimonsen@GMail.com
% Version 5.5
% 31/03-2015

% =========================================================================
% ===== Default values ====================================================
% =========================================================================

Rows = -1;
Cols = -1;
Scr = 'C';
Pos = 'C';
Tight = false;
Menu = '';
Transpose = false;
Extra = 'restart';
Offset = 0;
Not = [];
Only = [];
Simulink = 'include';
Adjust = [0,0];
Scale = 1;
Scale_Orig = [];
FreezeSize = false;
Order = [];  

% =========================================================================
% ===== Get inputs ========================================================
% =========================================================================

if (mod(nargin,2) ~= 0)
	error('Uneven number of inputs! Returning...');
end

for i = 1:2:nargin
	if (ischar(varargin{i+1}))
		varargin{i+1} = lower(varargin{i+1});
	end
	switch lower(varargin{i})
		case {'scr','screen','s'}
			% =============================================================
			Allow = {
				'C'		'c'		'center'	'main'		'primary'
				'W'		'w'		'west'		'l'			'left'
				'E'		'e'		'east'		'r'			'right'
				'N'		'n'		'north'		't'			'top'
				'S'		's'		'south'		'b'			'bottom'
				'Ext'	'ext'	'external'	'secondary'	NaN
				};
			Trigger = any(strcmp(Allow,varargin{i+1}),2);
			if (~any(Trigger))
				warning('Input ''%s'' not recognized for ''%s''!',varargin{i+1},varargin{i});
			else
				Scr = Allow{Trigger,1};
			end
		case {'pos','position','p'}
			% =============================================================
			Allow = {
				'C'		'c'		'center'	NaN		NaN
				'W'		'w'		'west'		'l'		'left'
				'E'		'e'		'east'		'r'		'right'
				'N'		'n'		'north'		't'		'top'
				'S'		's'		'south'		'b'		'bottom'
				'NW'	'nw'	NaN			NaN		NaN
				'NE'	'ne'	NaN			NaN		NaN
				'SW'	'sw'	NaN			NaN		NaN
				'SE'	'se'	NaN			NaN		NaN
				};
			Trigger = any(strcmp(Allow,varargin{i+1}),2);
			if (~any(Trigger))
				warning('Input ''%s'' not recognized for ''%s''!',varargin{i+1},varargin{i});
			else
				Pos = Allow{Trigger,1};
			end
		case {'rows','r'}
			% =============================================================
			if (strcmp(varargin{i+1},'auto'))
				Rows = -1;
			else
				if ~((isnumeric(varargin{i+1})) && (mod(varargin{i+1},1) == 0))
					warning('Input to ''%s'' must be a round numeric value!',varargin{i});
				else
					Rows = varargin{i+1};
				end
			end
		case {'cols','columns','c'}
			% =============================================================
			if (strcmp(varargin{i+1},'auto'))
				Cols = -1;
			else
				if ~((isnumeric(varargin{i+1})) && (mod(varargin{i+1},1) == 0))
					warning('Input to ''%s'' must be a round numeric value!',varargin{i});
				else
					Cols = varargin{i+1};
				end
			end
		case {'ext','extra'}
			% =============================================================
			Allow = {
				'ignore'
				'restart'
				};
			Trigger = any(strcmp(Allow,varargin{i+1}),2);
			if (~any(Trigger))
				warning('Input ''%s'' not recognized for ''%s''!',varargin{i+1},varargin{i});
			else
				Extra = Allow{Trigger,1};
			end
		case {'menu','bar'}
			% =============================================================
			Allow = {
				'figure'
				'none'
				};
			Trigger = any(strcmp(Allow,varargin{i+1}),2);
			if (~any(Trigger))
				warning('Input ''%s'' not recognized for ''%s''!',varargin{i+1},varargin{i});
			else
				Menu = Allow{Trigger,1};
			end
		case {'transpose'}
			% =============================================================
			if (~islogical(varargin{i+1}))
				warning('Input to ''%s'' must be logical!',varargin{i});
			else
				Transpose = varargin{i+1};
			end
		case {'offset','o'}
			% =============================================================
			if ~((isnumeric(varargin{i+1})) && (mod(varargin{i+1},1) == 0))
				warning('Input to ''%s'' must be a round numeric value!',varargin{i});
			else
				Offset = varargin{i+1};
			end
		case {'adjust','a'}
			% =============================================================
			if ~(isnumeric(varargin{i+1}))
				warning('Input to ''%s'' must be numeric!',varargin{i});
			elseif (numel(varargin{i+1}) ~= 2)
				warning('Input to ''%s'' must be a (1 x 2) matrix with integers or scalars less than one!',varargin{i});
			elseif (all(mod(varargin{i+1},1) == 0))
				Adjust = varargin{i+1};
			else
				if (any(varargin{i+1} > 1))
					warning('Input to ''%s'' must be either integers or scalars less than one!',varargin{i});
				else
					Adjust = varargin{i+1};
				end
			end
		case {'not'}
			% =============================================================
			if ~((isnumeric(varargin{i+1}) && all(mod(varargin{i+1},1) == 0)) || isobject(varargin{i+1}))
				warning('Input to ''%s'' must be round numeric values!',varargin{i});
			else
				if (isobject(varargin{i+1}))
					Not = get(varargin{i+1},'Number');
					if (iscell(Not))
						Not = cell2mat(Not);
					end
				else
					Not = varargin{i+1};
				end
			end
		case {'only'}
			% =============================================================
			if ~((isnumeric(varargin{i+1}) && all(mod(varargin{i+1},1) == 0)) || isobject(varargin{i+1}))
				warning('Input to ''%s'' must be round numeric values!',varargin{i});
			else
				if (isobject(varargin{i+1}))
					Only = get(varargin{i+1},'Number');
					if (iscell(Only))
						Only = cell2mat(Only);
					end
				else
					Only = varargin{i+1};
				end
			end
		case {'simu','simulink'}
			% =============================================================
			Allow = {
				'include'	NaN
				'exclude'	'ignore'
				'only'		NaN
				};
			Trigger = any(strcmp(Allow,varargin{i+1}),2);
			if (~any(Trigger))
				warning('Input ''%s'' not recognized for ''%s''!',varargin{i+1},varargin{i});
			else
				Simulink = Allow{Trigger,1};
			end
		case {'tight'}
			% =============================================================
			if (~islogical(varargin{i+1}))
				warning('Input to ''%s'' must be logical!',varargin{i});
			else
				Tight = varargin{i+1};
			end
		case {'scale'}
			% =============================================================
			if (~isnumeric(varargin{i+1}) || (numel(varargin{i+1}) > 2))
				warning('Input to ''%s'' must be one or two numeric values (x and y scale)!',varargin{i});
			else
				Scale = varargin{i+1};
			end
		case {'origin','orig','scaleorig','scaleorigin'}
			% =============================================================
			if (~isnumeric(varargin{i+1}) || (numel(varargin{i+1}) ~= 2))
				warning('Input to ''%s'' must be two numeric values (x0 and y0)!',varargin{i});
			else
				Scale_Orig = varargin{i+1};
			end
		case {'freezesize','freeze'}
			% =============================================================
			if (~islogical(varargin{i+1}))
				warning('Input to ''%s'' must be logical!',varargin{i});
			else
				FreezeSize = varargin{i+1};
			end
		case {'order'}
			% =============================================================
			if ~(isnumeric(varargin{i+1}) && all(mod(varargin{i+1},1) == 0))
				warning('Input to ''%s'' must be round numeric values!',varargin{i});
			else
				Order = varargin{i+1}; 
			end
		otherwise
			% =============================================================
			fprintf('Input ''%s'' not recognized!',varargin{i});
	end
end

% =========================================================================
% ===== Constants =========================================================
% =========================================================================

% ===== Limits ============================================================
Limits = {...
	'C'		[0,1]		[0,1]
	'W'		[0,0.5]		[0,1]
	'E'		[0.5,1]		[0,1]
	'N'		[0,1]		[0.5,1]
	'NW'	[0,0.5]		[0.5,1]
	'NE'	[0.5,1]		[0.5,1]
	'S'		[0,1]		[0,0.5]
	'SW'	[0,0.5]		[0,0.5]
	'SE'	[0.5,1]		[0,0.5]
	};

% ===== Screen angles =====================================================
Screens = {...
	(-pi * 1 / 4)	(+pi * 1 / 4)	'E'
	(-pi * 3 / 4)	(-pi * 1 / 4)	'S'
	(+pi * 1 / 4)	(+pi * 3 / 4)	'N'
	(+pi * 3 / 4)	Inf				'W'
	-Inf			(-pi * 3 / 4)	'W'
	};

% ===== Taskbar height ====================================================
Taskbar_Height = 40;

% =========================================================================
% ===== Figure list =======================================================
% =========================================================================

% ===== All figures =======================================================
Fig_Object = findall(0,'type','figure');
if (isempty(Fig_Object))
	fprintf('distFig: No figures to distribute!\n');
	return;
end

% ===== Remove not visible figures ========================================
Visible = arrayfun(@(n) (strcmp(get(Fig_Object(n),'Visible'),'on')),1:numel(Fig_Object));
Fig_Object = Fig_Object(Visible);

% ===== Include simulink figures ==========================================
Fig_Number = {Fig_Object.Number};
Empty = arrayfun(@(n) (isempty(Fig_Number{n})),1:numel(Fig_Number));
Fig_Number(Empty) = {0};
Fig_Number = cell2mat(Fig_Number);
[~,Index] = sort(get(Fig_Object(Fig_Number == 0),'Name'));
Temp = find(Fig_Number == 0);
Fig_Number(Temp(Index)) = -(numel(Index):(-1):1);

% ===== Logical distribution array ========================================
Fig_Dist = true(1,numel(Fig_Object));

% ===== Simulink figures ==================================================
if strcmp(Simulink,'exclude')
	Fig_Dist = (Fig_Number > 0);
end
if strcmp(Simulink,'only')
	Fig_Dist = (Fig_Number < 0);
end

% ===== Not ===============================================================
Fig_Dist(ismember(Fig_Number,Not)) = false;

% ===== Only ==============================================================
if (~isempty(Only))
	Fig_Dist(~ismember(Fig_Number,Only)) = false;
end

% ===== Number of figures =================================================
nFig_Dist = sum(Fig_Dist);

% ===== Return ============================================================
if (nFig_Dist == 0)
	return;
end

% =========================================================================
% ===== Monitor ===========================================================
% =========================================================================

% ===== Get monitor positions =============================================
Monitors = get(0,'MonitorPositions');
Monitors_Center = Monitors(:,1:2) + Monitors(:,3:4) / 2;
Monitors_Label = cell(size(Monitors,1),1);

% ===== Center monitor ====================================================
Temp = sum(abs(Monitors(:,1:2) - 1),2);
Trigger_C = find(Temp == min(Temp));
Monitors_Label{Trigger_C} = 'C';

% ===== External monitors =================================================
Dir = Monitors_Center - repmat(Monitors_Center(Trigger_C,:),size(Monitors_Center,1),1);
Theta = atan2(Dir(:,2),Dir(:,1));
for i = 1:numel(Theta)
	if (i ~= Trigger_C)
		Monitors_Label{i} = Screens{(Theta(i) > [Screens{:,1}]) & (Theta(i) <= [Screens{:,2}]),3};
	end
end

% ===== Select monitor ====================================================
Monitor = Monitors(strcmp(Monitors_Label,Scr),:);
if (isempty(Monitor))
	if (strcmp(Scr,'Ext'))
		if (size(Monitors,1) == 1)
			warning('No external monitor could be found!');
			Monitor = Monitors(1,:);
		else
			Monitor = Monitors(find(~strcmp(Monitors_Label,'C'),1,'first'),:);
		end
	else
		warning(['Screen ''%s'' could not be found! (Using ''C'' instead)\nThe following monitors were identified: %s.\n'...
			'If the function cannot find the monitor try and restart MATLAB - this will reset the monitor-matrix, which is being used by the function.'],...
			Scr,sprintf('%s\b\b\b\b\b',sprintf('''%s'' and ',Monitors_Label{:})));
		Monitor = Monitors(strcmp(Monitors_Label,'C'),:);
		if (isempty(Monitor))
			warning(['distFig cannot find matching monitors - using the first one!',...
				'\nPlease report this problem to AndersSSimonsen@GMail.com, where you include this matrix:'],0);
			disp(Monitors);
			Monitor = Monitors(1,:);
		end
	end
end

% =========================================================================
% ===== Rows and columns ==================================================
% =========================================================================

Limit = Limits(strcmp(Limits(:,1),Pos),:);
if ((Rows == (-1)) && (Cols == (-1)))
	AR_Ideal = 1.2;
	AR_Mean = Inf;
	for i = 1:50
		% ===== Exception =================================================
		if ((nFig_Dist + Offset) == 4)
			Rows = 2;
			Cols = 2;
			break;
		end
		
		% ===== Temporary rows and columns ================================
		Rows = i;
		Cols = ceil((nFig_Dist + Offset) / Rows);
		
		% ===== Sizes =====================================================
		x = round(linspace(Limit{2}(1),Limit{2}(2),Cols + 1) * Monitor(3)) + 1;
		y = round(linspace(Limit{3}(1),Limit{3}(2),Rows + 1) * (Monitor(4) - Taskbar_Height)) + Taskbar_Height + 1;
		clear Size;
		[Size(:,:,1),Size(:,:,2)] = meshgrid(diff(x),diff(y));
		
		% ===== Output ====================================================
		AR_Mean_Prev = AR_Mean;
		AR_Mean = mean(mean(Size(:,:,1) ./ Size(:,:,2)));
		if (AR_Mean > AR_Ideal)
			AR_Error = abs([AR_Ideal - AR_Mean_Prev,AR_Ideal - AR_Mean]);
			if (AR_Error(1) < AR_Error(2))
				Rows = i - 1;
			end
			Cols = ceil((nFig_Dist + Offset) / Rows);
			clear Size;
			break;
		end
	end
elseif (Rows == (-1))
	Rows = ceil((nFig_Dist + Offset) / Cols);
elseif (Cols == (-1))
	Cols = ceil((nFig_Dist + Offset) / Rows);
end

% =========================================================================
% ===== Location matrix ===================================================
% =========================================================================

Sets = ceil(nFig_Dist / (Rows * Cols) + Offset);
Mat = nan(Rows,Cols,Sets);
if (Transpose)
	Mat = permute(Mat,[2,1,3]);
end
if (isempty(Order))
	Order = (1:sum(Fig_Dist));
end
Temp = sort(Fig_Number(Fig_Dist));
if (max(Order) > numel(Temp))
	error('The ordering command is related to the input - not the actual figure numbers! Returning...');
end
Mat((1:nFig_Dist) + Offset) = Temp(Order);
if (Transpose)
	Mat = permute(Mat,[2,1,3]);
end
if (strcmp(Extra,'ignore'))
	Sets = 1;
end

% =========================================================================
% ===== Position and size =================================================
% =========================================================================

% ===== Grid ==============================================================
x = round(linspace(Limit{2}(1),Limit{2}(2),Cols + 1) * Monitor(3)) + 1;
y = round(linspace(Limit{3}(1),Limit{3}(2),Rows + 1) * (Monitor(4) - Taskbar_Height)) + Taskbar_Height + 1;

% ===== Scale =============================================================
if (isempty(Scale_Orig))
	Scale_Orig = [x(1),y(end)];
end
if ((all(mod(Scale_Orig,1) == 0)) && (any(Scale_Orig > 1)))
	ScaleOrig_Pixel = Scale_Orig;
else
	ScaleOrig_Pixel = Monitor(3:4) .* Scale_Orig;
end
if (numel(Scale) == 1)
	Scale(2) = Scale(1);
end
x = (x - ScaleOrig_Pixel(1)) * Scale(1) + ScaleOrig_Pixel(1); 
y = (y - ScaleOrig_Pixel(2)) * Scale(2) + ScaleOrig_Pixel(2); 

% ===== Positon ===========================================================
[FPos(:,:,1),FPos(:,:,2)] = meshgrid(x(1:end-1),y(1:end-1));
FPos(:,:,2) = flipud(FPos(:,:,2));
FPos(:,:,1) = FPos(:,:,1) + Monitor(1) - 1;
FPos(:,:,2) = FPos(:,:,2) + Monitor(2) - 1;

% ===== Size ==============================================================
[Size(:,:,1),Size(:,:,2)] = meshgrid(diff(x),diff(y));

% ===== Adjust ============================================================
if (any(Adjust ~= 0))
	if (all(mod(Adjust,1) == 0))
		Adjust_Pixels = Adjust;
	else
		Adjust_Pixels = round(Adjust .* (Monitor(3:4) + [0,-Taskbar_Height]));
	end
else
	Adjust_Pixels = [0,0];
end

% =========================================================================
% ===== Distribute figures ================================================
% =========================================================================

for s = 1:Sets
	for i = 1:size(Mat,1)
		for j = 1:size(Mat,2)
			if (~isnan(Mat(i,j,s)))
				% ===== Position and figure ===============================
				Pos_Temp = [squeeze(FPos(i,j,:))' + Adjust_Pixels,squeeze(Size(i,j,:))'];
				gcf_Temp = Fig_Object(Fig_Number == Mat(i,j,s));
				gca_Temp = get(gcf_Temp,'CurrentAxes');
				
				% ===== Freeze size =======================================
				if (FreezeSize)
					Pos_Orig = get(gcf_Temp,'Position');
					Pos_Temp(2) = Pos_Temp(2) + Pos_Temp(4) - Pos_Orig(4);
					Pos_Temp(3:4) = Pos_Orig(3:4);
				end
				
				% ===== Focus figures to bring to front ===================
				figure(gcf_Temp);
				
				% ===== Set units =========================================
				if (~isempty(Menu))
					set(gcf_Temp,'MenuBar',Menu);
				end
				if (Mat(i,j,s) < 0)
					drawnow;
				end
				gca_Units_Prev = get(gca_Temp,'Units');
				gcf_Units_Prev = get(gcf_Temp,'Units');
				
				% ===== Distribute figure =================================
				set(gca_Temp,'Units','normalized');
				set(gcf_Temp,'Units','pixels','OuterPosition',Pos_Temp);
				set(gcf_Temp,'Units',gcf_Units_Prev);
				
				% ===== Apply tight =======================================
				if ((Tight) && (Mat(i,j,s) > 0))
					drawnow;
					set(gca_Temp,'Position',[get(gca_Temp,'TightInset') * eye(4,2),1 - get(gca_Temp,'TightInset') * [1,0,1,0;0,1,0,1]']);
				end
				set(gca_Temp,'Units',gca_Units_Prev);
			end
		end
	end
end

% =========================================================================
% ===== Check for new updates =============================================
% =========================================================================

% This feature check for a new version of distFig. It is only an attempt,
% so it might not work without bugs. It can be disables by setting the
% variable below called "Update" to "false" or by deleting this whole
% section - this won't affect the rest of the function.
%
% It generates a .txt file called "distFig_Update.txt" located in the same
% folder as distFig.m, which holds the date for the last time the function
% checked for updates. In order to not check for updates all the time, the
% function only does this ONCE every 7 days (variable "Check_Interval"),
% and at max once for every started MATLAB session. Also, it uses appdata,
% so if you're using this as well, it might display messages multiple times
% within a session, if an update is available, if appdata is deleted.
%
% The update feature works by scanning the file exchange page for distFig
% for a string value called "datePublished", and compares it to the date of
% the current version (variable "Date_Current"). If these are not the same,
% a hyperlink will be shown to the function page on file exchange page. It
% thus depends on a very specific format, which might differ(?) from
% different languages - I'm not aware whether this is true, so please
% contact me on AndersSSimonsen@GMail.com if you find this to be true, so I
% can fix the feature. :)

% ===== DISABLE THE UPDATE FEATURE BY SETTING THIS VALUE TO false =========
Update = true;
% ===== DISABLE THE UPDATE FEATURE BY SETTING THIS VALUE TO false =========

if (Update)
	try
		% ===== Initialize ================================================
		Date_Current = [2015,03,31];
		Check_Interval = 7;
		Check = false;
		
		% ===== Path of update file =======================================
		Path = mfilename('fullpath');
		Path = sprintf('%s\\distFig_Update.txt',Path(1:end-numel('\distFig')));
		
		if (isempty(getappdata(0,'distFig_Update')))
			% ===== Set appdata ===========================================
			setappdata(0,'distFig_Update',1);
			
			% ===== Read update file ======================================
			File_ID = fopen(Path,'r');
			Clock = clock;
			Date_Today = Clock(1:3);
			if (File_ID == (-1))
				% ===== Generate update file ==============================
				File_ID = fopen(Path,'w');
				fprintf(File_ID,'%d %d %d',Date_Today(1),Date_Today(2),Date_Today(3));
				fclose(File_ID);
				Check = true;
			else
				% ===== Read last date ====================================
				Data = textscan(File_ID,'%s','delimiter','');
				if (strcmp(Data{1},'DISABLE'))
					Check = false;
				else
					Date_Last_Checked = cell2mat(textscan(Data{1}{1},'%f %f %f'));
					
					% ===== Write new date ================================
					File_ID = fopen(Path,'w');
					fprintf(File_ID,'%d %d %d',Date_Today(1),Date_Today(2),Date_Today(3));
					fclose(File_ID);
					
					% ===== Check for new updates? ========================
					Days_Elapsed = etime([Date_Today,zeros(1,3)],[Date_Last_Checked,zeros(1,3)]) / (60^2 * 24);
					if (Days_Elapsed > Check_Interval)
						Check = true;
					end
				end
			end
			
			% ===== Check for updates =====================================
			if (Check)
				% ===== Read File Exchange date ===========================
				Data = urlread('http://se.mathworks.com/matlabcentral/fileexchange/37176-distribute-figures');
				Temp = Data((-1:100) + regexp(Data,'datePublished'));
				Date_Online = cell2mat(textscan(Temp(26:end),'%f-%f-%f'));
				
				% ===== New update available? =============================
				if (~all(Date_Current == Date_Online))
					setappdata(0,'distFig_Update',2);
					setappdata(0,'distFig_Online_Date',Date_Online);
				end
			end
		end
		
		% ===== Message ===================================================
		if (getappdata(0,'distFig_Update') == 2)
			try
				% ===== New version =======================================
				Date_Online = getappdata(0,'distFig_Online_Date');
				fprintf('New update available for distFig! ');
				disp('<a href = "http://se.mathworks.com/matlabcentral/fileexchange/37176-distribute-figures">Download here</a>')
				fprintf('\tCurrent version date:\t%02.0f-%02.0f-%02.0f\n',Date_Current(1),Date_Current(2),Date_Current(3));
				fprintf('\tOnline version date:\t%02.0f-%02.0f-%02.0f\n',Date_Online(1),Date_Online(2),Date_Online(3));
				
				% ===== Postpone ==========================================
				fprintf('Click here to postpone this message for %d days: ',Check_Interval);
				Command = sprintf([...
					'setappdata(0,''distFig_Update'',1);',...
					'fprintf(''distFig update message postponed for %d days!\\n'');'
					],Check_Interval);
				disp(sprintf('<a href="matlab:%s">Postpone update message</a>',Command)); %#ok<DSPS>
				
				% ===== Disable ===========================================
				fprintf('Click here to disable the update feature forever: ');
				Command = sprintf([...
					'File_ID = fopen(''%s'',''w'');',...
					'fprintf(File_ID,''DISABLE'');',...
					'fclose(File_ID);',...
					'setappdata(0,''distFig_Update'',1);',...
					'fprintf(''distFig update feature disabled!\\n'');'
					],Path);
				disp(sprintf('<a href="matlab:%s">Disable update feature</a>',Command)); %#ok<DSPS>
			catch
				setappdata(0,'distFig_Update',1);
			end
		end
	catch Error %#ok<NASGU>
		% ===== Error handling ============================================
		try
			Path_Error = [Path(1:end-numel('\distFig_Update.txt')),'\distFig_Error.mat'];
			fprintf([
				'There was an error in distFig when checking for updates.',...
				'\nPlease report this to AndersSSimonsen@GMail.com by sending the error-file "distFig_Error.mat" located in %s. :)\n\n'],Path_Error);
			AppData = getappdata(0); %#ok<NASGU>
			save(Path_Error);
		catch
		end
	end
end