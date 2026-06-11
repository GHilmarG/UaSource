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
% as otherwise I'm just repeating the same output interval.
%
% In other parts of the code I use 
% 
%         ReminderFraction(CtrlVar.time,CtrlVar.DefineOutputsDt) < (CtrlVar.dt/(10*CtrlVar.DefineOutputsDt))) 
%
% to determine if I'm close enough to an output interval
% So if the remaining time interval from time to LimitTime is smaller than this, them I'm 
% at the same output interval.
%
%
% (dtOut/(10*CtrlVar.DefineOutputsDt)) ;
%
% If LimitTime is only fractionally larger than time, I may already have created an output file for this interval
% Currenlty, I'm not keeping track of the output intervals (maybe put this into RunInfo in the future?) So I'm having to 
% make a conservative guess.  If the difference is a small fraction of dtOut, then I will have created an output file
% already.  
% 
if  (LimitTime-time)<dtOut/1000
    LimitTime=LimitTime+Dt;
end


% Is it also possible that in absolut terms, dt is already on input very small
% so then also jump to next output interval

if abs(LimitTime-time) < 1000*eps
    LimitTime=LimitTime+Dt;
end

if (time+dtOut)> LimitTime  % if needed, redefine time step so that the LimitTime is not overstepped
    dtOut=LimitTime-time;   % (always strickly positive)
end

% Finally, check if the remaining time interval towards LimitTime is 
% a small fraction of dtOut. If so, then extend the time step all the way to LimitTime

RemainingDt=LimitTime-(time+dtOut);

if RemainingDt>0 && RemainingDt/dtOut < 1e-2
    dtOut=dtOut+RemainingDt;
end

end
