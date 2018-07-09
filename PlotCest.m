%load SlipPertEstimate
load SlipPertEstimateLooksFineToMeFor_m3.mat
H0=1000;
%%


vidObj = VideoWriter('CestLooksFineToMeFor_m3.avi'); open(vidObj);

figure('Position',[100 50 900 600],'units','pixel') ; 

Cfig=trisurf(TRIxy,x/H0,y/H0,Cest,'EdgeColor','none')  ; 
lightangle(-45,30) ; lighting phong ; shading interp
title('Estimated basal slipperinesss') ; 
xlabel('xps (km)') ; ylabel('yps (km)')  ; zlabel('(m kPa^{-3} a^{-1})')
axis tight
colorbar

for I=1:360
    
    view(I-37.5,30)
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    
end


close(vidObj);

% u=C tau^m , m/a=C kPa^m  => [C]=m kPa^m/a