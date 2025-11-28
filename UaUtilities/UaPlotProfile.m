

function UaPlotProfile(CtrlVar,MUA,F)


xyProfile=[-433.352739726027         -1024.47864625302 ; ...
            151.872078968574         -1491.31829170024 ];

xyProfile=xyProfile*1000;

x1=xyProfile(1,1); x2=xyProfile(2,1);
y1=xyProfile(1,2); y2=xyProfile(2,2);
nPoints=100;

xProfile=linspace(x1,x2,nPoints);  yProfile=linspace(y1,y2,nPoints);

Fs=scatteredInterpolant(F.x,F.y,F.s);
Fb=scatteredInterpolant(F.x,F.y,F.b);
FB=scatteredInterpolant(F.x,F.y,F.B);
Frho=scatteredInterpolant(F.x,F.y,F.rho);

sProfile=Fs(xProfile,yProfile);
bProfile=Fb(xProfile,yProfile);
BProfile=FB(xProfile,yProfile);

rhoProfile=Frho(xProfile,yProfile); 

hProfile=sProfile-bProfile;
SProfile=hProfile*0; 

% Due to interpolation, the floating conditions may not hold perfectly for the interpolated variables. 
% So, optionally, it might be best to recalculate s and b based on flotation along the profile. 

[bProfile,sProfile,hProfile]=Calc_bs_From_hBS(CtrlVar,MUA,hProfile,SProfile,BProfile,rhoProfile,F.rhow);


ProfileOffset=0; % sqrt( (xgl-xProfile(1)).^2 + (ygl-yProfile(1)).^2 ) ;
Profile=sqrt( (xProfile-xProfile(1)).^2 + (yProfile-yProfile(1)).^2 )-ProfileOffset;


fig=FindOrCreateFigure("UaPlotProfile") ; clf(fig)


plot(Profile/1000,BProfile,'k'); hold on
plot(Profile/1000,bProfile,'b');
plot(Profile/1000,sProfile,'b');

% Bedrock polygon
BxPoly=[min(Profile)/1000  max(Profile)/1000  fliplr(Profile)/1000 ] ;
ByPoly=[-2000  -2000   fliplr(BProfile) ] ;
fill(BxPoly,ByPoly,[128 128 128]/255) ;

OceanxPoly=[Profile/1000 fliplr(Profile)/1000 ] ;
OceanyPoly=[BProfile   fliplr(bProfile) ] ;
fill(OceanxPoly,OceanyPoly,[0 0 1]) ;


ICExPoly=[Profile/1000 fliplr(Profile)/1000 ] ;
ICEyPoly=[bProfile   fliplr(sProfile) ] ;
fill(ICExPoly,ICEyPoly,[0.58 0.815 0.988]) ;
title("Profile",Interpreter="latex",FontSize=16)


axis tight
ylim([min(BProfile) max(sProfile)]) ;

xlabel("distance (km)",Interpreter="latex") ; ylabel("(m)",Interpreter="latex")




end