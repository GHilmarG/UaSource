
function ID=MaskAntarctica(xloc,yloc)

persistent mask x y


if isempty(mask)


    filename="NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc";


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