
function ID=MaskAntarctica(xloc,yloc,filename)

%% For each location over Antarctica, returns:
%
%
%  0 = ocean
%  1 = ice-free land
%  2 = grounded
%  3 = floating ice
%  4 = Lake Vostok
%
% xloc and yloc can be vectors
%
% This is based on finding the closest (x,y) location in the mask array to the (xloc,yloc). 
%
% There is no interpolation used, therefore the returned values are always the integer values of the mask.
%
% The mask here is that used by Bedmachine, but it would be easy to adapt this for any other mask.
%
% filename : file with mask, for example "NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc"
%
%%


persistent mask x y


if isempty(mask)

    if nargin < 3 || isempty(filename)
        filename="NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc";
    end


    %% Read data

    fprintf(" Reading BedMachine data from file %s ",filename)
    x = ncread(filename,'x');
    y = ncread(filename,'y');
    % bed = ncread(filename,'bed')';
    % surface = ncread(filename,'surface')';
    % thickness = ncread(filename,'thickness')';
    % firn = ncread(filename,'firn')';
    mask = ncread(filename,'mask')'; mask=double(mask) ;
    % source = ncread(filename,'source')';
    % geoid = ncread(filename,'geoid')'; geoid=double(geoid) ;
    % errbed = ncread(filename,'errbed')';
    y=flipud(y);
    mask=rot90(mask,3);
    % figure ; imagesc(x,y,mask) ; axis equal xy tight ; colorbar ; clim([-0.5 4.5]) ; colormap(jet(5))
end

dx=double(x(2)-x(1)) ; xmin=double(min(x)) ; 
dy=double(y(2)-y(1)) ; ymin=double(min(y)) ; 


ix=round((xloc-xmin)/dx) +1;
iy=round((yloc-ymin)/dy) +1;

ind=sub2ind(size(mask),ix,iy);

ID=mask(ind);

%  ocean/ice/land mask
%  0 = ocean
%  1 = ice-free land
%  2 = grounded
%  3 = floating ice
%  4 = Lake Vostok
%
end