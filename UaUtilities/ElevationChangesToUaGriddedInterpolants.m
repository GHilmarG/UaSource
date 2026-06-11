
%%
%
% Schröder, Ludwig; Horwath, Martin; Dietrich, Reinhard; Helm, Veit; van den Broeke, Michiel R; Ligtenberg, Stefan R M (2019): 
% Gridded surface elevation changes from multi-mission satellite altimetry 1978-2017. PANGAEA, 
% https://doi.org/10.1594/PANGAEA.897390,

load('MeshBoundaryCoordinatesForAntarcticaBasedOnBedmachine','Boundary')
load('GroundingLineForAntarcticaBasedOnBedmachine','GroundingLines') ; 

filename="sec_mrg.nc" ; 
finfo=ncinfo(filename);
fprintf(" Reading data from file... %s ",filename) 
x = ncread(filename,'x');           x=double(x);
y = ncread(filename,'y');           y=double(y);
time = ncread(filename,'time');     time=double(time) ; 
sec = ncread(filename,'sec');
sec_std = ncread(filename,'sec_std');
fprintf("...done.\n")

% I'm confused about what is actually being shown. I think it is the absolute elevation changes
% in meters since september 2010, which appears to be about in sec(:,:,386)

%% Some plotting to check

% figure ; contour(x,y,sec(:,:,end)) ;

close all

figure
for I=386-100:10:474  % 1:50:474
    
    contourf(x/1000,y/1000,sec(:,:,I)',LineStyle="none") ; 
    axis equal  ; cbar=colorbar ; caxis([-50 5])
    title(sprintf("elevation change between %f and Sept 2010",time(I)))
    xlabel("xps (km)")
    ylabel("yps (km)")
    title(cbar,"(m)")
    pause(1)

    %prompt = 'Do you want more? Y/N [Y]: ';
    %str = input(prompt,'s');
    %if isempty(str)
    %    str = 'Y';
    %end
end

%%
% Create xyt interpolant
Fdh=griddedInterpolant({x,y,time},sec) ;


%% Some further plots
Year=1985:10:2015 ;

Year=2000;
dh=Fdh({x,y,Year});

figure
contourf(x/1000,y/1000,sec(:,:,I)',LineStyle="none") ;
axis equal  ; cbar=colorbar ; caxis([-50 5])
title(sprintf("elevation change between %f and Sept 2010",Year))
xlabel("xps (km)")
ylabel("yps (km)")
title(cbar,"(m)")

%% Quick estimate of dh/dt
Year1=2010 ;
Year2=2018 ;

dh=Fdh({x,y,[Year1; Year2]});

dhdt=(dh(:,:,2)-dh(:,:,1))/(Year2-Year1);

figure
contourf(x/1000,y/1000,dhdt',LineStyle="none") ;
axis equal  ; cbar=colorbar ; caxis([-8 1])
title(sprintf("rate of elevation change between %f and %2",Year1,Year2))
xlabel("xps (km)")
ylabel("yps (km)")
title(cbar,"(m/yr)")
%%

xloc=-1580e3 ; yloc=-230e3;  % approx PIG
[~,ix]=min(abs(x-xloc));
[~,iy]=min(abs(y-yloc));
Y1=1990 ; Y2=2020 ; 
I=find(time>=Y1 & time <=Y2) ;

figure ;
t=time(I); dhloc=sec(ix,iy,I) ; dhloc=dhloc(:);
I=~isnan(dhloc) ; t=t(I) ; dhloc=dhloc(I);
plot(t,dhloc,'o-r')
ylabel('Elevation change (m)')
xlabel('time (years)')
title(sprintf("Elevation change at PS=(%i,%i) (km) \n with respect to sep 2010",xloc/1000,yloc/1000))

N=numel(dhloc) ; A=[ones(N,1) t]; b=dhloc ; sol=A\b;

hold on
plot(t,sol(1)+sol(2)*t,'k')
text((Y1+0.1*(Y2-Y1)),0,sprintf("dh/dt=%f (m/yr)",sol(2)))


%% 
ixRange=10:250; iyRange=50:350;

ixRange=1:numel(x); iyRange=1:numel(y);

Y1=2000 ; Y2=2018 ;
iTimeRange=find(time>=Y1 & time <=Y2) ;
t=time(iTimeRange); t=t(:) ;
dhdt=nan(numel(x),numel(y)) ;

for ix=ixRange
    %fprintf("%i\n",ix)
    for iy=iyRange

        
        dhloc=sec(ix,iy,iTimeRange) ; 
        dhloc=dhloc(:);

        I=~isnan(dhloc) ;   % get rid of NaNs 
        tt=t(I) ; 
        dhloc=dhloc(I); 

        N=numel(dhloc) ;
        if N>5
            A=[ones(N,1) tt];
            b=dhloc ; sol=A\b;
            dhdt(ix,iy)=sol(2);
        end
    end
end
[xGL,yGL]=ReadBindschadlerGroundingLine;
%%
fig=FindOrCreateFigure("dh/dt");
Nlevels=250; 
contourf(x/1000,y/1000,dhdt',Nlevels,LineStyle="none") ;
axis equal  ; cbar=colorbar ; caxis([-5 1])
title(sprintf("Mean rate of elevation change between %i and %i \n (Schröder et al 2019: https://doi.org/10.1594/PANGAEA.897390),",Y1,Y2))
xlabel("xps (km)")
ylabel("yps (km)")
title(cbar,"(m/yr)")
hold on
plot(xGL/1000,yGL/1000,'k')
plot(GroundingLines(:,1)/1000,GroundingLines(:,2)/1000,'k')
colormap(othercolor('RdYlBu_11b',2*Nlevels))
ModifyColormap
hold on ; PlotLatLonGrid(1000,5,45,2000,"black",true);
axis off


%%  Now combine this with estimates of dh/dt over the ice shelves
% First step is to get a mask for floating areas
% Here use the mask from Bedmachine
%  ocean/ice/land mask
%  0 = ocean
%  1 = ice-free land
%  2 = grounded
%  3 = floating ice 
%  4 = Lake Vostok
filename = "../Bedmachine/BedMachineAntarctica_2020-07-15_v02.nc" ;
fprintf(" Reading BedMachine data from file %s ",filename) 
xBM = ncread(filename,'x'); xBM=double(xBM);
yBM = ncread(filename,'y');  yBM=double(yBM);
%bedBM = ncread(filename,'bed')'; 
%surfaceBM = ncread(filename,'surface')'; 
%thicknessBM = ncread(filename,'thickness')'; 
%firnBM = ncread(filename,'firn')'; 
maskBM = ncread(filename,'mask')'; maskBM=double(maskBM) ;
%source = ncread(filename,'source')'; 
%geoid = ncread(filename,'geoid')'; geoid=double(geoid) ; 
%errbed = ncread(filename,'errbed')'; 
fprintf("...done.\n")




FBmask=griddedInterpolant({double(xBM),double(flipud(yBM))},rot90(double(maskBM),3)) ;
Mask=FBmask({x,y});

FindOrCreateFigure("Bedmachine mask")
levels=5 ; 
contourf(x/1000,y/1000,Mask',levels,LineStyle="none") ;
colorbar ; caxis([0 4]) ; axis equal
colormap jet

%%  % now read in dh/dt over ice shelves
%
% 
% Hi Hilmar,
% 
% Here is a summary of the various files:
% 
% https://www.dropbox.com/s/0r0nshdaqavdyqh/ice_shelf_ra_v2.h5?dl=0
% 
% 
% 
% This file contains 1994-2017 data from ERS-1, ERS-2, Envisat, and CryoSat-2 in regions where data was available from all four radar altimeters. Data is up to 81.5S only. Thickness change is in m/y averaged over 1994-2017.
% 
% https://www.dropbox.com/s/t9bnx48ve0s5dph/ice_shelf_ra_v2.1-81.5S.h5?dl=0
% 
% This file contains average thickness change (in m/y) averaged over 1994-2017 using ERS-1, ERS-2, Envisat, and CryoSat-2 data upto latitude 81.5S. Below that latitude, data is average thickness change (also in m/y) from 2010-2017 (CryoSat-2 data only). 
% 
% https://www.dropbox.com/s/rr5nvh6l56glwqz/ice_shelf_ra_v2.1-cs2-only.h5?dl=0 
% 
% 
% This file has complete Antarctic ice shelf coverage from 2010-2017 (only CryoSat-2 data was used). Thickness change is in m/y averaged over 2010-2017.
% 
% Best,
% Susheel
%
% reading data
Version="v2" ;   % 1994-2017 data, only up to 81.5 S
Version="v2.1"; % don't use
Version="v2.1-81.5S"; % Same as v2 but with CryoSat-2 south of 81.5S, has a jumb and is questionable to use
Version="v2.1-cs2-only"; % complete 2010-2017 coverage using only CryoSat-2, average thickness change in m/y averaged over 2010-2017
filename="ice_shelf_ra_"+Version+".h5";  % 
filename="../IceShelfThinningRates/Suadusum/"+filename; 
h5disp(filename)

% consider using either 1) all date from 1994-2017 and only north of 81.5 and 2)
% all data from 2010 to 2017 which is only CryoSat-2,
% focus on 1994-2017 compilation (only north of 81.5)

lat = h5read(filename,'/lat');
lon = h5read(filename,'/lon');
thickness_change = h5read(filename,'/thickness_change');

%time = h5read(filename,'/time');


% creating data vectors in xps yps space
[Lat,Lon]=ndgrid(lat,lon);
[xps,yps]=ll2xy(Lat(:),Lon(:));

temp=thickness_change' ; temp=temp(:) ; thickness_change=temp ;

% simple visualisation to see if this makes sense
FindOrCreateFigure("Thickness change, point data")
plot3(xps/1000,yps/1000,thickness_change,'.') ; title(' thickness change ' )

FdhdtIceShelves=scatteredInterpolant(xps,yps,thickness_change) ;

dhdtIceShelves=FdhdtIceShelves({x,y});

FindOrCreateFigure("dh/dt over ice shelves")
levels=20 ; 
contourf(x/1000,y/1000,dhdtIceShelves',levels,LineStyle="none") ;
colorbar ; caxis([-4 1]) ; axis equal

%%

dhdtAll=dhdt ;
I=isnan(dhdtIceShelves) ; dhdtIceShelves(I)=0 ; 

% Note: There is so little data on dh/dt over the iceshelves for the Amundsen Sea
% that after inspection, it appears questionable to try to combine it
% I here simply put some resonable made up dh/dt values for the Amundsen sea ice shelves
%
PIGIceShelfdhdt=-2 ;
TWGIceShelfdhdt=-2 ; 
SKIceShelfdhdt=-3 ; 
[X,Y]=ndgrid(x,y) ;

I=X>-1700e3 & X < -1500e3 &  Y>-380e3 & Y < -200e3 ;   dhdtIceShelves(I)=PIGIceShelfdhdt ;  %PIG
I=X>-1700e3 & X < -1500e3 &  Y>-550e3 & Y < -400e3 ;   dhdtIceShelves(I)=TWGIceShelfdhdt ;  % TWG
I=X>-1700e3 & X < -1300e3 &  Y>-900e3 & Y < -550e3 ;   dhdtIceShelves(I)=SKIceShelfdhdt ;  % Smith and rest

I=Mask==3; dhdtAll(I)=dhdtIceShelves(I);   % floating ice shelves 
I=Mask==0; dhdtAll(I)=0;   %  ocean
I=Mask==1; dhdtAll(I)=0;   %  ice free land
I=isnan(dhdtAll); dhdtAll(I)=0;   %  get rid of NaNs


fig=FindOrCreateFigure("dh/dt combined");
Nlevels=250; 
contourf(x/1000,y/1000,dhdtAll',Nlevels,LineStyle="none") ;
axis equal  ;
title(sprintf("Mean rate of elevation change between %i and %i \n (Schröder et al 2019: https://doi.org/10.1594/PANGAEA.897390) \n Susheel et al,",Y1,Y2))
xlabel("xps (km)")
ylabel("yps (km)")

hold on
%plot(xGL/1000,yGL/1000,'k')


hold on ; PlotLatLonGrid(1000,5,45,2000,"black",true);
hold on ; plot(Boundary(:,1)/1000,Boundary(:,2)/1000,'k')
plot(GroundingLines(:,1)/1000,GroundingLines(:,2)/1000,'k')
axis off
colormap(othercolor('RdYlBu_11b',Nlevels))
cbar=colorbar ; caxis([-5 1]) ; title(cbar,"(m/yr)")
ModifyColormap




%% Checking

load MUA_Antarctica.mat; 
xFE=MUA.coordinates(:,1) ; yFE=MUA.coordinates(:,2) ; 

dhdt=Fdh2000to2018(xFE,yFE);

FindOrCreateFigure("dh/dt over mesh")
PlotMeshScalarVariable([],MUA,dhdt);
cbar=colorbar ; caxis([-5 1]) ; title(cbar,"(m/yr)")
ModifyColormap


%%
% Save

Fdh2000to2018=griddedInterpolant({x,y},dhdtAll) ;

save("FdhdtMeasuredRatesOfElevationChanges2000to2018","Fdh2000to2018")


