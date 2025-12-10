





function [Fs,Fb,FB,Frho]=BedMachineGreenlandToUaGriddedInterpolants(filename,N,SaveOutputs) 

%%  
% Reads BedMachine nc file and creates Ua gridded interpolants with s, b, B and rho.
% 
% Note: This is just a rough idea as to how this could be done and most likely will need some modifications for some particular
% applications.
%
%
% filename              name nc of file with the BedMachine data
%
% N                     subsample options, use every N data point in x and y direction. (N=1 implies no sub-sampling.)
%
%
%
% Examples:   
%
% Using default values: 
%
%       BedMachineToUaGriddedInterpolants ; 
%
%  
% Maximum data resolution (N=1) and boundary resolution equal to the resolution of the data sets (500 meters). 
%   
%       BedMachineToUaGriddedInterpolants("BedMachineAntarctica_2020-07-15_v02.nc",1,false,500,true);
%
% Reduced resolution (N=10), boundary resolution at 20km, and do not save outputs.
%
%       BedMachineToUaGriddedInterpolants("BedMachineAntarctica_2020-07-15_v02.nc",10,false,20e3,false) ; 
%
% Full data resolution (N=1), boundary resolution at 1000m and saving outputs.
%   
%       BedMachineToUaGriddedInterpolants("BedMachineAntarctica_2020-07-15_v02.nc",1,false,1000,true);
%
%
% NOTE: The BEDMACHINE data compilation assumes ocean density of 1027 kg/m^3. So you MUST set 
%       rhow=1027 in your own user input DefineDensity.m or DefineGeometryAndDensities.m m-files.
%
%% 



arguments
    filename (1,1) string = "BedMachineGreenland-v5.nc";
    N (1,1) double = 2
    SaveOutputs (1,1) logical = true ; 
end


%% Read data

fprintf(" Reading BedMachine data from file %s ",filename) 



finfo=ncinfo(filename);
x = double(ncread(filename,'x'));
y = double(ncread(filename,'y'));
bed = ncread(filename,'bed'); 
surface = ncread(filename,'surface'); 
thickness = ncread(filename,'thickness'); 
mask = ncread(filename,'mask');  mask=double(mask) ;



source = ncread(filename,'source'); 
geoid = ncread(filename,'geoid'); geoid=double(geoid) ; 
errbed = ncread(filename,'errbed'); 
fprintf("...done.\n")
% All heights are referenced to mean sea level using the geod EIGEN-6D4
%
% To obtain height with respect to WGS84 add geod
%
%
%% Plot raw data, just to see if all is OK
%
% fprintf(' Generating some figures of raw data...')
% figure(10) ; imagesc(x,y,bed); axis xy equal; caxis([-4000 1800]); title(' bed ' ) ; colorbar ; axis tight
% figure(20) ; imagesc(x,y,thickness); axis xy equal; caxis([0 4500]); title(' thickness ' ) ; colorbar ; axis tight
% figure(40) ; imagesc(x,y,surface); axis xy equal; caxis([0 4800]); title(' surface ' ) ; colorbar ; axis tight
% figure(50) ; imagesc(x,y,geoid); axis xy equal; caxis([-100 100]); title(' geoid ' ) ; colorbar ; axis tight
% figure(60) ; imagesc(x,y,mask); axis xy equal; caxis([0 4]); title(' mask ' ) ; colorbar ; axis tight
% drawnow ; fprintf('done.\n')


%% Check relationship between geometrical variables
% figure(100) ; imagesc(x,y,(mask~=0).*((surface-(bed+thickness+firn)))); axis xy equal;  title('(mask~=0).*(surface-(bed+thickness+firn))' ) ; colorbar ; axis tight ; caxis([-50 50])
% figure(110) ; imagesc(x,y,(mask~=0).*(surface-(bed+thickness))); axis xy equal;  title('(mask~=0).*(surface-(bed+thickness))' ) ; colorbar ; axis tight ; caxis([-50 50])
% figure(120) ; imagesc(x,y,(firn+thickness)); axis xy equal;  title('firn+thickness' ) ; colorbar ; axis tight ;

%% Ua variables
%
% 
s=surface ;
b=surface-thickness ;
h=thickness ; % or just h=s-b 
B=bed ; 



OceanMask=(mask==0) ; 
CanadaMask=(mask==4) ;

rho=s*0+917 ;% is this the value they used?!
rhow=1027;



%%
%% Possible subsampling

x=x(1:N:end) ;
y=y(1:N:end) ; 
s=s(1:N:end,1:N:end) ;
rho=rho(1:N:end,1:N:end) ;
b=b(1:N:end,1:N:end) ;
B=B(1:N:end,1:N:end) ;
mask=mask(1:N:end,1:N:end) ;
[X,Y]=ndgrid(x,y) ;


%%


fprintf(' Plotting s, b, h, B and rho over the data grid')
figure(200) ; imagesc(x,y,s'); axis xy equal; clim([0 4000]); title(' s ' ) ; colorbar ; axis tight
figure(210) ; imagesc(x,y,b'); axis xy equal; clim([-2000 4000]); title(' b ' ) ; colorbar ; axis tight
figure(220) ; imagesc(x,y,h'); axis xy equal; clim([0 4000]); title(' h ' ) ; colorbar ; axis tight
figure(230) ; imagesc(x,y,B'); axis xy equal; clim([-4000 4000]); title(' B ' ) ; colorbar ; axis tight
figure(240) ; imagesc(x,y,rho'); axis xy equal; clim([0 920]); title(' rho ' ) ; colorbar ; axis tight
drawnow ; fprintf('done.\n')


%% Boundaries


fBoundary=FindOrCreateFigure("mask boundaries") ; clf(fBoundary)
mask2=mask;
mask2(mask2==4)=0;   % getting rid of the Canadian land mass that just causes problems, in mask2 this is now an ocean!  
[M]=contour(X,Y,mask2,[0.5 1.5 2.5,3.5]); 
colorbar
axis equal


% Coast lines
level=0.5 ;  % this contour level must be in M
I=find(M(1,:)==level) ; [~,J]=max(M(2,I)) ; 
[Nglpoints,iloc]=sort(M(2,I),'descend');
fprintf(' %i points in the longest contour line segment.\n',M(2,I(J)) );
Ngl=20 ;   % a bit random, but taking just these longest line segments provides sufficient information for plotting purposes. 
for k=1:Ngl
    GLvector{k}=M(:,I(iloc(k))+1:I(iloc(k))+M(2,I(iloc(k)))) ;
end

GL=[]; 
for k=1:Ngl
    temp=GLvector{k}';
    GL=[GL ; temp ; nan nan];
end
CoastLine=GL;



% Grounded Ice boundary
level=1.5 ;  % this contour level must be in M
I=find(M(1,:)==level) ; [~,J]=max(M(2,I)) ; 
[Nglpoints,iloc]=sort(M(2,I),'descend');
fprintf(' %i points in the longest contour line segment.\n',M(2,I(J)) );
Ngl=20 ;    % a bit random, but taking just these longest line segments provides sufficient information for plotting purposes. 
for k=1:Ngl
    GLvector{k}=M(:,I(iloc(k))+1:I(iloc(k))+M(2,I(iloc(k)))) ;
end

GL=[]; 
for k=1:Ngl
    temp=GLvector{k}';
    GL=[GL ; temp ; nan nan];
end
IceBoundary=GL;


fCO=FindOrCreateFigure("Boundaries") ; clf(fCO)
plot(IceBoundary(:,1)/1000,IceBoundary(:,2)/1000,'b')  ; 
hold on
plot(CoastLine(:,1)/1000,CoastLine(:,2)/1000,'k')  ; 
axis equal ;
axis([-645 857 -3370 -645])
legend("Ice boundary","Coastlines")

% save("GreenlandBoundaries.mat","CoastLine","IceBoundary")


%% create gridded interpolants

fprintf(' Creating gridded interpolants...')
[X,Y]=ndgrid(double(x),double(flipud(y))) ;

Fs=griddedInterpolant(X,Y,rot90(double(s'),3)); 
Fb=griddedInterpolant(X,Y,rot90(double(b'),3)); 
FB=griddedInterpolant(X,Y,rot90(double(B'),3)); 
Frho=griddedInterpolant(X,Y,rot90(double(rho'),3)); 
Fmask=griddedInterpolant(X,Y,rot90(double(mask'),3)); 

fprintf('done.\n')

if SaveOutputs
    GriddedFileName=replace(filename+"-Stride"+num2str(N)+"-GriddedInterpolantes.mat",".nc","");
    fprintf(' Saving BedMachineGreenlandGriddedInterpolants with Fs, Fb, FB and Frho in %s. \n',GriddedFileName)
    save(GriddedFileName,"Fs","Fb","FB","Frho","Fmask","CoastLine","IceBoundary","-v7.3")
end

%% Test gridded interpolant

% Here creating a simple square mesh
fprintf(' Testing interpolants by mapping on a FE grid...')
CtrlVar=Ua2D_DefaultParameters();
CtrlVar.UaSquareMesh.xmin=min(x(:));
CtrlVar.UaSquareMesh.xmax=max(x(:));
CtrlVar.UaSquareMesh.ymin=min(y(:));
CtrlVar.UaSquareMesh.ymax=max(y(:));

CtrlVar.UaSquareMesh.nx=100;
CtrlVar.UaSquareMesh.ny=round((CtrlVar.UaSquareMesh.ymax-CtrlVar.UaSquareMesh.ymin)/(CtrlVar.UaSquareMesh.xmax-CtrlVar.UaSquareMesh.xmin))*CtrlVar.UaSquareMesh.nx;

[~,~,MUA]=UaSquareMesh(CtrlVar);
FindOrCreateFigure("MUAMESH") ; PlotMuaMesh([],MUA)


%load("MUA_Greenland.mat","MUA")


CtrlVar.PlotXYscale=1000;
F=UaFields;
F.x=MUA.coordinates(:,1) ; 
F.y=MUA.coordinates(:,2) ; 



F.s=Fs(F.x,F.y);
F.b=Fb(F.x,F.y);
F.B=FB(F.x,F.y);
F.rho=Frho(F.x,F.y);
F.rhow=rhow;
F.S=F.s*0; 

Mask=Fmask(F.x,F.y) ; 
[F.b,F.h,F.GF]=Calc_bh_From_sBS(CtrlVar,MUA,F.s,F.B,F.S,F.rho,F.rhow) ; 



[X,Y]=meshgrid(x,y) ; [lat,lon]=psn2ll(X,Y);  % I need to use meshgrid here because I use contour for lat/lon



cbar=UaPlots(CtrlVar,MUA,F,F.s,FigureTitle="s") ;
colormap(cmocean('ice',25))  ;
xlabel('(km)' ) ; xlabel('(km)' ) ; title('b') ;
title(cbar,'m a.s.l')
hold on ;
plot(CoastLine(:,1)/1000,CoastLine(:,2)/1000,'w')  ;
plot(IceBoundary(:,1)/1000,IceBoundary(:,2)/1000,'b')  ;
LatLonGrid(X/1000,Y/1000,lat,lon,LineColor=[0.5 0.5 0.5],LabelSpacing=200);

cbar=UaPlots(CtrlVar,MUA,F,F.b,FigureTitle="b") ;
xlabel('(km)' ) ; xlabel('(km)' ) ; title('b') ;
title(cbar,'m a.s.l')
hold on ;
plot(CoastLine(:,1)/1000,CoastLine(:,2)/1000,'w')  ;
plot(IceBoundary(:,1)/1000,IceBoundary(:,2)/1000,'b')  ;
LatLonGrid(X/1000,Y/1000,lat,lon,LineColor=[0.5 0.5 0.5],LabelSpacing=200);
CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);

cbar=UaPlots(CtrlVar,MUA,F,F.B,FigureTitle="B") ;
xlabel('(km)' ) ; xlabel('(km)' ) ; title('B') ;
title(cbar,'m a.s.l')
hold on ; LatLonGrid(X/1000,Y/1000,lat,lon,LineColor=[0.5 0.5 0.5],LabelSpacing=200);
colormap(othercolor("Mdarkterrain",25))  ;
clim([-1500 2800])


cbar=UaPlots(CtrlVar,MUA,F,Mask,FigureTitle="Mask") ;
xlabel('(km)' ) ; xlabel('(km)' ) ; title('Mask') ;
title(cbar,'m a.s.l')
hold on ; LatLonGrid(X/1000,Y/1000,lat,lon,LineColor=[0.5 0.5 0.5],LabelSpacing=200);
clim([-0.5 4.5])
colormap(othercolor("Mdarkterrain",5))  ;


fprintf('done.\n')

return
end







