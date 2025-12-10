


%%

H=1; eta=1/2; rho=917. ;

sigma_cx=10 ; sigma_cy=10 ; ampl_c=0;
sigma_sx=10 ; sigma_sy=10 ; ampl_s=0;
theta=0.*pi/180.;


t=NaN;


% [-------Experiment F
alpha=3*pi/180. ; ampl_b=0.1 ; ampl_c=0. ; sigma_bx=10 ; sigma_by=10 ; C=1 ;
nx=2^10 ; ny=nx ; dx= 400/nx ; dy=dx ;
cmin=-0.02 ; cmax=0.04 ; levels=[cmin:(cmax-cmin)/20:cmax];
%-----------]


% % [------- Small Gauss, high slip
% alpha=3*pi/180. ; ampl_b=0.1 ; ampl_c=0. ; sigma_bx=0.5 ; sigma_by=0.5 ; C=100 ;
% nx=1024 ; ny=nx ; dx= 100/nx ; dy=dx ;
% cmin=-0.015 ; cmax=0.1 ; levels=[cmin:(cmax-cmin)/30:cmax];
% %-----------]


g=1/(rho*sin(alpha)*H);


x=[-dy*ny/2:dy:dy*(ny/2-1)]';  % model domain
y=[-dy*ny/2:dy:dy*(ny/2-1)]';  % model domain



% C  : mean slipperiness
% ampl_b and ampl_c : amplitudes of b and c perturbations
% sigma_bx          : width of Gaussian-shaped b perturbation, etc.
%

[s,u,v,w,db,dc]=Gauss_xy(x,dx,y,dy,alpha,C,ampl_b,sigma_bx,sigma_by,ampl_c,sigma_cx,sigma_cy,theta);


% plot surface topography

%cmin=min(s(:)); cmax=max(s(:));
%levels=[cmin:(cmax-cmin)/30:cmax];



axis equal ; axis tight
[C1,h1]=contourf(x,y,s',levels);
colorbar
xlabel('x (h)') ; ylabel('y (h)') ; title('GS Surface topography')
caxis([cmin cmax])


figure

speed=sqrt(u.*u+v.*v);
levels2=[0:0.01:0.1];
[C4,h4]=contour(x,y,speed',levels2);
caxis([0 0.1])
set(h4,'ShowText','on','TextStep',0.0)
hold on
istep=8;
h5=quiver(x(1:istep:end),y(1:istep:end),u(1:istep:end,1:istep:end)',v(1:istep:end,1:istep:end)');
axis equal tight
axis([-30 30 -30 30])
set(h5,'Color','k')
xlabel('x (h)') ; ylabel('y (h)') ; title('Flow over a Gaussian peak');
%colorbar
%axis equal




