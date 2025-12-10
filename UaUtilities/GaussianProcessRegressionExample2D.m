

%%
%
% Example of how to use Gaussian Process Regression model on two-dimentional (spatial) data
%




nx=50;  ny=nx; 
x=linspace(-5,5,nx); 
y=linspace(-5,5,ny); 
[X,Y]=ndgrid(x,x) ;

Z=exp(-X.*X-Y.*Y) + 0.01*randn(nx,ny)  ;




xy=[X(:) Y(:)] ; z=Z(:); 

tic
%gprMdl=fitrgp(xy,z,PredictorNames={'x','y'},ResponseName="z");
gprMdl=fitrgam(xy,z,PredictorNames={'x','y'},ResponseName="z",FitStandardDeviation=true,OptimizeHyperparameters="auto");
toc

% Here the prediction is evaluated at the data locations
tic
[zpred,~,zint] = predict(gprMdl,xy);
toc

Zpred=reshape(zpred,nx,ny);




figure(100)
surf(x,y,Z)
title("data")


figure(150)
surf(x,y,Zpred)
title("fit to data")

Residuals=Zpred-Z; 

figure(200) ; surf(x,y,Residuals)  ; title("Residuals")


figure(300) ; histogram(Residuals)


%%