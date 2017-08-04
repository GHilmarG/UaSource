 

F=ExplicitEstimationForUaFields(CtrlVar,F,F0,Fm1);


[F1.ub,F1.vb,F1.ud,F1.vd,F1.h]=...
    ExplicitEstimation(CtrlVar.dt,CtrlVar.dtRatio,CtrlVar.CurrentRunStepNumber,...
    F0.ub,F0.dubdt,Fm1.dubdt,...
    F0.vb,F0.dvbdt,Fm1.dvbdt,...
    F0.ud,F0.duddt,Fm1.duddt,...
    F0.vd,F0.dvddt,Fm1.dvddt,...
    F0.h,F0.dhdt,Fm1.dhdt);




%%
clear all
fVector=[1  ;2  ;3   ;4  ;  5 ; 6 ] ;

t=      [1  ; 1.2 ; 3 ; 5.09999 ;  5.1 ; 5.7  ; 7 ; 8.5 ; 9.3] ;
t=      [1  ; 1.2 ; 3 ; 5.0 ;  5.1 ; 5.10000001  ; 5.1002 ; 8.5 ; 9.3] ;
fVector=t.^2; dfdtExact=2*t;
%fVector=sin(t*2*pi/10/2);

%fVector=t;

dfdtm1=[]; dfdt=[]; dtm1=[]; dtRatio=[];

Itime=6;

dt=t(Itime+1)-t(Itime);
f0=fVector(Itime) ;
if Itime>1
    dfdt=(fVector(Itime)-fVector(Itime-1))./(t(Itime)-t(Itime-1));
    dtm1=t(Itime)-t(Itime-1);
    f1m=fVector(Itime-1);
end

if Itime>2
    dfdtm1=(fVector(Itime-1)-fVector(Itime-2))./(t(Itime-1)-t(Itime-2));
end

if Itime>1
    dtRatio=dt/dtm1;
end

dfdt=dfdtExact(Itime);  dfdtm1=dfdtExact(Itime-1); 


extrapolation='constant';
extrapolation='linear';
extrapolation='quadratic';

f1=ExplicitEstimation2(dt,dtRatio,Itime,extrapolation,f0,dfdt,dfdtm1);

fprintf('\n\n f1=%g \t %g\n',f1,fVector(Itime+1))