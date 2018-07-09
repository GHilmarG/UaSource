
load SSS2dforward Experiment s h b B a u v w etaInt time connectivity coordinates Nele nod nip Nnodes AGlen C ...
    Luv Luvrhs lambdau Lh Lhrhs lambdah n m alpha rhow rho g

sMeas=real(s);
uMeas=real(u);
vMeas=real(v);
wMeas=real(w);
bMeas=real(b);
BMeas=real(B);
xMeas=coordinates(:,1);
yMeas=coordinates(:,2);

%save SyntDataCPert sMeas uMeas vMeas wMeas bMeas BMeas xMeas yMeas
save SyntDataBPert sMeas uMeas vMeas wMeas bMeas BMeas xMeas yMeas