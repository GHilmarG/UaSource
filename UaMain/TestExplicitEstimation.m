
%%
%
% Simple test of the explicit estimation used in Ua
% 
%  dy/dt=f(t)
%
% here:
%
%  f(t)=sin(t) with f(t=0)=0; 
%
% Note: this is a particularly simple test, and more generally we expect
%
% dy/dt=f(y(t),t)
%
% But this is a sufficient test to make sure that no typos were introduced in the implementation of the AB2 for variable time
% steps.
%
%%

ExplicitEstimatorDegree=2; % 0 : use last value as an estimate for next
                           % 1 : linear extrapolation (i.e. forward Euler)
                           % 2 : second order Adams-Bashforth extrapolation (for a variable time-step)

% Running this code for ExplicitEstimatorDegree=1  and then again with ExplicitEstimatorDegree=2, clearly demonstrates the
% increased accuracy obtained with the two-step second-order Adams-Bashforth expression.

N=200;  % Number of time steps


y=nan(N,1);


y(1)=0 ;   % Initial condition 

dy0dt=nan;
T=4*2*pi; 
t=linspace(0,T,N)';                  % constant time step
t=cumsum(rand(N,1)) ; t=T*(t-t(1))/max(t) ; % random time step

%t(1:3)=t(1:3)/10; % Initially the method is linear, and there is an accumulation of error in the first few time steps

dt=nan; 

for Itime=1:N-1

    % t(Itime) is the "current time", ie the time for which the solution is already known/estimated.
    %
    % An estimate for y at time t(Itime+1) is being sought. The new estimate is stored at y(Itime+1)
    % 
    % As this is an explicit estimate f(t)=sin(t) is only evaluated at the current time, i.e. at time t(Itime)
    %
    %


    dym1dt=dy0dt ; % previous dy/dt values, ie from time t(Itime-1)
    y0=y(Itime)  ; % y0 is the estimate of the solution at current time. This is the future estimate (y(Itime+1)) from previous iteration
    dy0dt=sin(t(Itime)) ;
    dtm1=dt;
    dt=t(Itime+1)-t(Itime);
    dtRatio=dt/dtm1;
    
    IDegree=min(Itime,ExplicitEstimatorDegree+1);
    y(Itime+1)=ExplicitEstimation(dt,dtRatio,IDegree,y0,dy0dt,dym1dt);


end


yAnalytical=-cos(t)+1; 

figure(10) ; 
yyaxis left
hold off ;  
plot(t,y,"b.",LineWidth=2,DisplayName="Numerical") ; 
hold on ; 
plot(t,yAnalytical,"k-",DisplayName="Analytical")
ylabel("Numerical and Analytical Solutions")

yyaxis right
Err=y-yAnalytical;
hold off
plot(t,Err,"r--",DisplayName="Numerical-Analytical")
ylabel("Error")


ld=legend;
xlabel("time")
axis tight

%%