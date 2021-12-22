function dtOut=NoOverStepping(CtrlVar,time,dtIn,Dt)
%
%
%  On input Dt is the maximum time step, while dtIn is the current time step
%
% To ensure the dtOut is not larger than Dt one could simpy set
%  
%   dtOut=min(dtIn,Dt)
%
%  However, I also would like dtOut to be some simple fraction on Dt
%
%  I do this by defining a LimitTime as 
%
%   LimitTime=Dt*ceil(time/Dt)
%
%  and then I ensure that on return time+dtOut does not exceed this limit time. 
%

%dtOut=round(dtIn,14); % round to 10 significant digits.
%                      % this is needed so that the sum time+dtOut is always numerically different from time
dtOut=dtIn ;


LimitTime=Dt*ceil(time/Dt);  % time that should not be overstepped.
                             % LimitTime is always >= time
                             
% If LimitTime is numerically effectivly equal to time, then advance LimitTime by Dt
if abs(LimitTime-time) < 1000*eps
    LimitTime=time+Dt;
end

if (time+dtOut)> LimitTime  % if needed, redefine time step so that the LimitTime is not overstepped
    dtOut=LimitTime-time;   % (always strickly positive)
end

% Finally, check if the remaining time interval towards LimitTime is 
% a small fraction of dtOut. If so, then extend the time step all the way to LimitTime

RemainingDt=LimitTime-(time+dtOut);

if RemainingDt/dtOut < 1e-2
    dtOut=dtOut+RemainingDt;
end

end
