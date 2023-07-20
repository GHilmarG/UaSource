function [Fs,Fb,FB,Frho,Boundary]=BedMachineToUaGriddedInterpolants(filename,N,LakeVostok,BoundaryResolution,SaveOutputs) 

%%  
% Reads BedMachine nc file and creates Úa gridded interpolants with s, b, B and rho.
% 
% Note: This is just a rough idea as to how this could be done and most likely will need some modifications for some particular
% applications.
%
%
% filename              name nc of file with the BedMachine data
%
% N                     subsample options, use every N data point in x and y direction. (N=1 implies no sub-sampling.)
%
% LakeVostok            logical variable, if set to false, then Lake Vostok is eliminated from the data by raising the bedrock to match the
%                       lower ice surface over the lake.
%
% BoundaryResolution    resolution in meters when contouring mask to produce boundary coordinates. This can not be smaller than the
%                       resolution of the data set itself, which is 500m. And if N is set to, for example N=10, then the smallest
%                       possible resolution for the boundary is N*500.
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
    filename (1,1) string = "BedMachineAntarctica_2020-07-15_v02.nc" ;
    N (1,1) double = 1
    LakeVostok (1,1) logical = false
    BoundaryResolution (1,1) double = 2*N*500;      % this mush be larger than 500m, and preferrably some interger multiple of that
    SaveOutputs (1,1) logical = false ; 
end


%% Read data

fprintf(" Reading BedMachine data from file %s ",filename) 
x = ncread(filename,'x');
y = ncread(filename,'y');
bed = ncread(filename,'bed')'; 
surface = ncread(filename,'surface')'; 
thickness = ncread(filename,'thickness')'; 
firn = ncread(filename,'firn')'; 
mask = ncread(filename,'mask')'; mask=double(mask) ;
source = ncread(filename,'source')'; 
geoid = ncread(filename,'geoid')'; geoid=double(geoid) ; 
errbed = ncread(filename,'errbed')'; 
fprintf("...done.\n")
% All heights are referenced to mean sea level using the geod EIGEN-6D4
%
% To obtain height with respect to WGS84 add geod
%
% Rema= Surface + Firn + Geod
%
%  ocean/ice/land mask
%  0 = ocean
%  1 = ice-free land
%  2 = grounded
%  3 = floating ice 
%  4 = Lake Vostok
%
%
% Apparantly  surface=bed+thickness over the grounded areas
%
% However this surface is not the REMA surface
%
% So the actual surface (with respect to sea level) is: s = surface+firn
%                                                       b=bed=surface-thickness=s-firn-thickness
%                                           therefore   h=s-b=surface+firn-(surface-thickness)=firn+thickness
%
%
%% Plot raw data, just to see if all is OK
%
fprintf(' Generating some figures of raw data...')
figure(10) ; imagesc(x,y,bed); axis xy equal; caxis([-4000 1800]); title(' bed ' ) ; colorbar ; axis tight
figure(20) ; imagesc(x,y,thickness); axis xy equal; caxis([0 4500]); title(' thickness ' ) ; colorbar ; axis tight
figure(30) ; imagesc(x,y,firn); axis xy equal; caxis([0 50]); title(' firn ' ) ; colorbar ; axis tight 
figure(40) ; imagesc(x,y,surface); axis xy equal; caxis([0 4800]); title(' surface ' ) ; colorbar ; axis tight
figure(50) ; imagesc(x,y,geoid); axis xy equal; caxis([-100 100]); title(' geoid ' ) ; colorbar ; axis tight
figure(60) ; imagesc(x,y,mask); axis xy equal; caxis([0 4]); title(' mask ' ) ; colorbar ; axis tight
drawnow ; fprintf('done.\n')


%% Check relationship between geometrical variables
% figure(100) ; imagesc(x,y,(mask~=0).*((surface-(bed+thickness+firn)))); axis xy equal;  title('(mask~=0).*(surface-(bed+thickness+firn))' ) ; colorbar ; axis tight ; caxis([-50 50])
% figure(110) ; imagesc(x,y,(mask~=0).*(surface-(bed+thickness))); axis xy equal;  title('(mask~=0).*(surface-(bed+thickness))' ) ; colorbar ; axis tight ; caxis([-50 50])
% figure(120) ; imagesc(x,y,(firn+thickness)); axis xy equal;  title('firn+thickness' ) ; colorbar ; axis tight ;

%% Ua variables
%
% I'm assuming here that the user wants the variable s to be surface+firn, in which case the vertically averaged ice density must be
% modfied accordingly.
%
% 
s=surface+firn ;
b=surface-thickness ;
h=firn+thickness ; % or just h=s-b 
B=bed ; 

if ~LakeVostok
    % get rid of Lake Vostok by setting lower ice surface (b) equal to bedrock elevation (B) over the lake as defined my mask.
    IV=mask== 4;
    B(IV)=b(IV);
end

OceanMask=(mask==0) ; 

rhoi=917 ; 
rhoMin=100; 


%%
%
% The 'firn' field in the Bedmachine data is not the firn thickness! It is the firn air content
%
% $$ \mathrm{firn} =  \frac{1}{\rho_i}  \int_{firn layer} (\rho_i - \rho_f) \, dz $$
%
% The 'firn air content' has the units distance.
%
% \rho = (1-firn/h) \rho_i   where h is the total ice thickness (Ua definition).
%
%
%%


rho=(1-firn./(h+eps)) .*rhoi ;  % vertically averaged ice density over the deph h=s-b=firn+thicknes

% Note: one would expect that the minimum value of rho would be equal to the firn density,
%       but there are a number of values where rho=0, that is when firn=h. Not sure how this can happen
%       when the ice thickness is larger than zero... Possibly some inconsistencies between the data sets providing air content (firn) and
%       ice thickness. Note: firn is 'air content', but 'air content'
I=h>0 ; figure(190) ; histogram(rho(I),'Normalization','probability') ; title('histogram of densities where h>0')

fprintf('Over glaciated areas %f%% of densities are smaller than %f kg/m^3 \n', 100*numel(find(rho(I)<rhoMin))/numel(find(I)),rhoMin)
fprintf('These densities are set to %f kg/m^3\n',rhoMin)
I=h>0 & rho<rhoMin ; rho(I)=rhoMin ; 

rho(OceanMask)=rhoi;  % be carefull here! To lessen the risk of potential extrapolation errors, I fill this with the ice dencities over the ocean.

fprintf(' Plotting s, b, h, B and rho over the data grid')
figure(200) ; imagesc(x,y,s); axis xy equal; caxis([0 4000]); title(' s ' ) ; colorbar ; axis tight
figure(210) ; imagesc(x,y,b); axis xy equal; caxis([-2000 4000]); title(' b ' ) ; colorbar ; axis tight
figure(220) ; imagesc(x,y,h); axis xy equal; caxis([0 4000]); title(' h ' ) ; colorbar ; axis tight
figure(230) ; imagesc(x,y,B); axis xy equal; caxis([-4000 4000]); title(' B ' ) ; colorbar ; axis tight
figure(240) ; imagesc(x,y,rho); axis xy equal; caxis([0 920]); title(' rho ' ) ; colorbar ; axis tight
drawnow ; fprintf('done.\n')
%% Possible subsampling

x=x(1:N:end) ;
y=y(1:N:end) ; 
s=s(1:N:end,1:N:end) ;
rho=rho(1:N:end,1:N:end) ;
b=b(1:N:end,1:N:end) ;
B=B(1:N:end,1:N:end) ;
mask=mask(1:N:end,1:N:end) ;

%% create gridded interpolants

fprintf(' Creating gridded interpolants...')
[X,Y]=ndgrid(double(x),double(flipud(y))) ;

Fs=griddedInterpolant(X,Y,rot90(double(s),3)); 
Fb=griddedInterpolant(X,Y,rot90(double(b),3)); 
FB=griddedInterpolant(X,Y,rot90(double(B),3)); 
Frho=griddedInterpolant(X,Y,rot90(double(rho),3)); 
fprintf('done.\n')

if SaveOutputs
    fprintf(' Saving BedMachineGriddedInterpolants with Fs, Fb, FB and Frho. \n')
    save('BedMachineGriddedInterpolants','Fs','Fb','FB','Frho','-v7.3')
end

%% Test gridded interpolants

% load in an old mesh of Antarctica and map on this mesh. Just done for testing purposes and to see if the gridded data looks sensible.
fprintf(' Testing interpolants by mapping on a FE grid of Antartica...')

load MUA_Antarctica.mat; 
xFE=MUA.coordinates(:,1) ; yFE=MUA.coordinates(:,2) ; 

sFE=Fs(xFE,yFE);
bFE=Fb(xFE,yFE);
BFE=FB(xFE,yFE);
rhoFE=Frho(xFE,yFE);

CtrlVar.PlotXYscale=1000;
 
bfig=FindOrCreateFigure('b') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,bFE);
xlabel('xps (km)' ) ; xlabel('yps (km)' ) ; title('b') ; title(cbar,'m a.s.l')


sfig=FindOrCreateFigure('s') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,sFE);
xlabel('xps (km)' ) ; xlabel('yps (km)' ) ; title('s') ; title(cbar,'m a.s.l')


Bfig=FindOrCreateFigure('B') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,BFE);
xlabel('xps (km)' ) ; xlabel('yps (km)' ) ; title('B') ; title(cbar,'m a.s.l')


rhofig=FindOrCreateFigure('rho') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,rhoFE); 
xlabel('xps (km)' ) ; xlabel('yps (km)' ) ; title('rho') ; title(cbar,'kg/m^3')

rhofig=FindOrCreateFigure('bFE-BFE') ; 
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,bFE-BFE); 
xlabel('xps (km)' ) ; xlabel('yps (km)' ) ; title('b-B') ; title(cbar,'m')

fprintf('done.\n')
%% Meshboundary coordinates
% Find the boundary by extracting the 0.5 contour line of the IceMask Then extract the largest single contourline. Depending on the
% situation, this may or may not be what the user wants. But this appears a reasonable guess as to what most users might want most of the
% time. 
fprintf('Creating MeshBoundary coordinates...')
IceMask=(mask~=0) ;  IceMask=double(IceMask);
DataResolution=N*500 ; % the resolution of the data set is 500 m
NN=min([round(BoundaryResolution/DataResolution) 1]) ;


fc=FindOrCreateFigure('contour') ;  
[M]=contour(x(1:NN:end),y(1:NN:end),IceMask(1:NN:end,1:NN:end),1) ; axis equal
hold on; plot(M(1,:), M(2, :), 'r.');

% now find longest contourline
level=0.5 ;  % this contour level must be in M
I=find(M(1,:)==level) ; [~,J]=max(M(2,I)) ; 
fprintf(' %i points in the longest contour line segment.\n',M(2,I(J)) );
Boundary=M(:,I(J)+1:I(J)+M(2,I(J))) ; 
plot(Boundary(1,:),Boundary(2,:),'o-k')  ; axis equal

Boundary=Boundary';
fprintf('done.\n')
if SaveOutputs
    
    fprintf('Saving MeshBoundaryCoordinates. \n ')
    save('MeshBoundaryCoordinatesForAntarcticaBasedOnBedmachine','Boundary')
end




%% Grounding line coordinates
% Find the boundary by extracting the 0.5 contour line of the IceMask Then extract the largest single contourline. Depending on the
% situation, this may or may not be what the user wants. But this appears a reasonable guess as to what most users might want most of the
% time. 
fprintf('Creating grounding line coordinates...')
GroundedIceMask=(mask==2 ) ;  GroundedIceMask=double(GroundedIceMask);

% N=10 ; BoundaryResolution=20e3;
DataResolution=N*500 ; % the resolution of the data set is 500 m
NN=min([round(BoundaryResolution/DataResolution) 1]) ;

%NN=1; % testing
fc=FindOrCreateFigure('GL contour') ;  
M=contour(x(1:NN:end),y(1:NN:end),GroundedIceMask(1:NN:end,1:NN:end),1) ; axis equal
hold on; plot(M(1,:), M(2, :), 'r.');

% now find Ngl longest contourline
level=0.5 ;  % this contour level must be in M
I=find(M(1,:)==level) ; [~,J]=max(M(2,I)) ; 

[Nglpoints,iloc]=sort(M(2,I),'descend');

fprintf(' %i points in the longest contour line segment.\n',M(2,I(J)) );

Ngl=2000 ;   % 2000 longest grounding lines
for k=1:Ngl
    GLvector{k}=M(:,I(iloc(k))+1:I(iloc(k))+M(2,I(iloc(k)))) ;
end

GL=[]; 
for k=1:Ngl
    temp=GLvector{k}';
    GL=[GL ; temp ; nan nan];
end

FindOrCreateFigure("Grounding lines")
plot(GL(:,1)/1000,GL(:,2)/1000,'r')  ; axis equal ;

fprintf('done.\n')

if SaveOutputs
    xGL=GL(:,1); yGL=GL(:,2);
    fprintf('Saving grounding line \n ')
    save('GroundingLineForAntarcticaBasedOnBedmachine','xGL','yGL') ; 
end




%%  Testing mesh boundary coordinates and creating  a new one with different spacing between points and some level of smoothing

CtrlVar.GLtension=1e-12; % tension of spline, 1: no smoothing; 0: straight line
CtrlVar.GLds=5e3 ; 


[xB,yB,nx,ny] = Smooth2dPos(Boundary(:,1),Boundary(:,2),CtrlVar);
MeshBoundaryCoordinates=[xB(:) yB(:)] ;
fc=FindOrCreateFigure('MeshBoundaryCoordinates') ;  
plot(MeshBoundaryCoordinates(:,1),MeshBoundaryCoordinates(:,2),'.-') ; axis equal
title("Example of a smoothed and resampled boundary")







