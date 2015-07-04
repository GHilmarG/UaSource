function dtOut=NoOverStepping(CtrlVar,time,dtIn,Dt)
%
%
%  Dt is a time interval and time+dtOut should not exceed this limit
%
%

dtOut=round(dtIn,14); % round to 10 significant digits.
                      % this is needed so that the sum time+dtOut is always numerically different from time

LimitTime=Dt*ceil(time/Dt);  % time that should not be overstepped.
                             % LimitTime is always >= time
                             
% If LimitTime is equal to time, then advance LimitTime by Dt
if abs(LimitTime-time) < CtrlVar.dtmin;
    LimitTime=time+Dt;
end

if (time+dtOut)> LimitTime  % if needed, redefine time step so that the LimitTime is not overstepped
    dtOut=LimitTime-time;   % (always strickly positive)
end

% Check if for the time step dtOut, the remaining time interval towards LimitTime is 
% a small fraction of this time step. If so, then extend the time step all the way to LimitTime

RemainingDt=LimitTime-(time+dtOut);

if RemainingDt/dtOut < 1e-2
    dtOut=dtOut+RemainingDt;
end

end
