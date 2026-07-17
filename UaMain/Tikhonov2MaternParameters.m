






function [alphaMatern,tauMatern,kappaMatern,sigma2Matern,nuMatern,rhoMatern]=Tikhonov2MaternParameters(ga,gs,Area)




alphaMatern=1;
tauMatern=gs/sqrt(Area);

if ga==0
    kappaMatern=0;
else
    kappaMatern=ga/(gs+eps(gs));
end

d=2; 
nuMatern=alphaMatern-d/2;

sigma2Matern = gamma(nuMatern) / (gamma(alphaMatern) * (4*pi)^(d/2) * kappaMatern^(2*nuMatern) * tauMatern^2);
rhoMatern=sqrt(8*nuMatern)/kappaMatern;


end