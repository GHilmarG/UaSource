function varargout=ExplicitEstimation(dt,dtRatio,Itime,varargin)

%%
%
% varargout=ExplicitEstimation(dt,dtRatio,Itime,varargin)
% explicit estimate using second-order variable time step Adams-Bashforth
% http://lucan.me.jhu.edu/wiki/index.php/Second-order_variable_time_step_Adams-Bashforth_method
%
%
% Two step Adams-Bashforth method, variable time step
%
% $$\frac{du}{dt}=f(u(t),t)$$
% 
%
% $$ u_{n+1}=u_n + \frac{\Delta t_n}{2 \Delta t_{n-1}} \left ( ( 2 \Delta t_{n-1} + \Delta t_n ) f(t_n,u_n) - \Delta t_n f(t_{n-1},u_{n-1}) \right ) $$
% 
% $$\Delta t_{n-1} = t_n-t_{n-1} $$
%
% $$\Delta t_{n} = t_{n+1}-t_{n} $$
%
% Can also be written as:
%
%
% $$ u_{n+1}=u_n + \Delta t_n \left ( f_n + \frac{1}{2} \frac{\Delta t_n}{\Delta t_{n-1}} \left ( f_n - f_{n-1} \right ) \right )  $$
% 
%
% and as:
%
% $$ u_{n+1}=u_n + \Delta t_n \left ( (1+ r/2)  f_n - r f_{n-1}/2 \right )   $$
%
% where
%
% $$ r:= \Delta t_n/\Delta t_{n-1} $$
%
% See Eq (17) in :
%
% Interval methods of Adams-Bashforth type with variable step sizes
% Andrzej Marciniak1,2 ·Malgorzata A. Jankowska3
% https://doi.org/10.1007/s11075-019-00774-y
%
% The constant time step expression when $\Delta t_{n}= \Delta t_{n-1}$, is 
%
% $$ u_{n+1}=u_n + \Delta t_n \left ( \frac{3}{2}   f(t_n,u_n) - \frac{1}{2} f(t_{n-1},u_{n-1}) \right ) $$
%
%
% varargin is in triples, where
%
% F0.h,F0.dhdt,Fm1.dhdt
%
% [h0, dh0/dt, dhm1/dt]
%
%
% See also: TestExplicitEstimation.m 
%
%%


nInputs=nargin-3;
nOutputs = nargout;
varargout = cell(1,nOutputs);

if nInputs~=3*nOutputs
    error(' wrong number of inputs ')
end

if Itime==1 

    for I=1:nOutputs
        varargout{I}=varargin{1+3*(I-1)};
    end
    
elseif Itime==2 
    
    for I=1:nOutputs
        varargout{I}=varargin{1+3*(I-1)}+dt*varargin{2+3*(I-1)};
    end


elseif Itime >= 3  % second-order variable time step Adams-Bashforth
   
 
    
    
    for I=1:nOutputs
        varargout{I}=varargin{1+3*(I-1)}+dt*((1+0.5*dtRatio)*varargin{2+3*(I-1)}-0.5*dtRatio*varargin{3+3*(I-1)});
    end


end

end