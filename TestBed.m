%%
%clear all ; close all 
x=linspace(0,1800e3);
y=linspace(-100e3,100e3);



[X,Y]=meshgrid(x,y);

sigma=5e3; 
sigma=10e3; 
w=30e3; 
d=1e3;

B=720-778.5*X/750e3+d*(1./(1+exp(-2*(Y-w)/sigma))+1./(1+exp(2*(Y+w)/sigma)));

figure
surfc(X,Y,B);

%% transverse profile

B=d*(1./(1+exp(-2*(y-w)/sigma))+1./(1+exp(2*(y+w)/sigma)));
%B=1./(1+exp(2*(y+w)/sigma));
figure ; plot(y,B)




% transverse slope

dBdy=-2*d*exp(-2*(y-w)/sigma)./(1+exp(-2*(y-w)/sigma)).^2/sigma+...
      2*d*exp(+2*(y+w)/sigma)./(1+exp(+2*(y+w)/sigma)).^2/sigma;
  
  figure ; plot(y,dBdy) ; title('transverse slope')
  
  



