


%%
Klear
c=0.5 ; Func=@(x)  (x-c).^2 ; slope0=-2*c; 
% c=0.75 ; Func=@(x)  (x-c).^2 ; slope0=-2*c; 
% c1=0.1 ; c2=0.75  ; Func=@(x)  (x-c1).^2 .* (x-c2).^2 ; slope0=-2*c1.* (0-c2).^2 - (0-c1).^2 .* 2.*c2;

% c1=50 ; c2=50  ; Func=@(x)  (x-c1).^2 .* (x-c2).^2 ; slope0=-2*c1.* (0-c2).^2 - (0-c1).^2 .* 2.*c2;

c1=0 ; c2=0.5  ; c3=1.1 ; Func=@(x)  -(x-c1).* (x-c2).*(x-c3) ; slope0=-(c1*c2+c1*c3+c2*c3);  f0=Func(0) ; f1=Func(1) ; 



[gmin,fmin,BackTrackInfo]=BackTracking(slope0,1,f0,f1,Func);

xvector=linspace(0,1) ; yvector=Func(xvector);

figure
plot(xvector,yvector) ; hold on
plot(gmin,fmin,'or')
plot(BackTrackInfo.InfoVector(:,1),BackTrackInfo.InfoVector(:,2),'+b')
%%





