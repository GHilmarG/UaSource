function [RunInfo,dtOut,dtRatio]=AdaptiveTimeStepping(UserVar,RunInfo,CtrlVar,MUA,F)

%% dtOut=AdaptiveTimeStepping(time,dtIn,nlInfo,CtrlVar)
%  modifies time step size
%
% Decision about increasing the time step size is based on the number of non-linear interations over the last few time steps.
%
% The main idea is to limit the number of non-linear iteration so that the NR is within the quadradic regime
% Experience has shown that a target number of iterations (CtrlVar.ATSTargetIterations) within 3 to 5 is good for this purpose
%
% Time step is increased if r<1 where
%
%     r=N/M
%
%   where
%   N is the max number of non-linear iteration over last n time steps
%   M is the target number of iterations
%
%   here
%     M=CtrlVar.ATSTargetIterations
%   and
%     n=CtrlVar.ATSintervalUp
%
%   (N does not need to be specifed.)
%
%  Currently the time step is only decreased if either:
%        a) the number of non-linear iterations in last time step was larger than 25
%        b) number of iterations over last n times steps were all larger than 10
%  where n=CtrlVar.ATSintervalDown
%
% There are some further modifications possible:
%  -time step is adjusted so that time interval for making transient plots (CtrlVar.TransientPlotDt) is not skipped over
%  -time step is not increased further than the target time step CtrlVar.ATSdtMax
%  -time step is adjusted so that total simulation time does not exceed CtrlVar.TotalTime
%
%
%

narginchk(5,5)

persistent dtNotUserAdjusted dtOutLast dtModifiedOutside

RunInfo.Forward.AdaptiveTimeSteppingTimeStepModifiedForOutputs=0 ;

time=CtrlVar.time;
dtIn=CtrlVar.dt ;



if isempty(dtModifiedOutside)
    dtModifiedOutside=false;
end




% potentially dt was previously adjusted for plotting/saving purposes
% if so then dtNotUserAdjusted is the previous unmodified time step
% and I do want to revert to that time step. I do so only
% however if dt has not been changed outside of function
if ~isempty(dtOutLast)
    if dtIn==dtOutLast  % dt not modified outside of function
        if ~isempty(dtNotUserAdjusted)
            dtIn=dtNotUserAdjusted ;
        end
    else

        RunInfo.Forward.AdaptiveTimeSteppingResetCounter=0;
        dtModifiedOutside=true;
    end
end

% dtOut=dtIn ;
dtOut=max(dtIn,CtrlVar.ATSdtMin) ;



% I first check if the previous forward calculation did not converge. If it did
% not converge I reduced the time step and reset all info about previous
% interations to reset the adaptive-time stepping approuch


RunInfo.Forward.AdaptiveTimeSteppingResetCounter=RunInfo.Forward.AdaptiveTimeSteppingResetCounter+1;


if CtrlVar.UaRunType=="-uvh-" % The time stepping algorithm is based on the number of uvh iterations, and hence only works for uvh

    if ~RunInfo.Forward.uvhConverged || dtModifiedOutside


        dtModifiedOutside=false ;

        if ~RunInfo.Forward.uvhConverged
            dtOut=dtIn/CtrlVar.ATStimeStepFactorDownNOuvhConvergence;
            fprintf(CtrlVar.fidlog,' ---------------- Adaptive Time Stepping: time step decreased from %-g to %-g due to lack of convergence in last uvh step. \n ',dtIn,dtOut);
        end


    elseif CtrlVar.AdaptiveTimeStepping && CtrlVar.CurrentRunStepNumber>1



        % I base the decision on the values in ItVector
        % Extract last

        nStepBacks=max([CtrlVar.ATSintervalDown,CtrlVar.ATSintervalUp]);
        ItVector=RunInfo.Forward.uvhIterations(max(CtrlVar.CurrentRunStepNumber-nStepBacks,1):CtrlVar.CurrentRunStepNumber-1);
        nItVector=numel(ItVector) ;


        % TimeStepUpRatio the ratio between maximum number of non-linear iterations over
        % last CtrlVar.ATSintervalUp iterations, divided by CtrlVar.ATSTargetIterations
        % It TimeStepUpRatio is smaller than 1, the number of non-linear iterations has
        % consistently been below target and time step should potentially be increased.

        if (CtrlVar.CurrentRunStepNumber-1)>=CtrlVar.ATSintervalUp
            TimeStepUpRatio=max(RunInfo.Forward.uvhIterations(max(CtrlVar.CurrentRunStepNumber-1-CtrlVar.ATSintervalUp+1,1):CtrlVar.CurrentRunStepNumber-1))/CtrlVar.ATSTargetIterations ;
        else
            TimeStepUpRatio=NaN ;
        end

        if (CtrlVar.CurrentRunStepNumber-1)>=CtrlVar.ATSintervalDown
            TimeStepDownRatio=min(RunInfo.Forward.uvhIterations(max(CtrlVar.CurrentRunStepNumber-1-CtrlVar.ATSintervalDown+1,1):CtrlVar.CurrentRunStepNumber-1))/CtrlVar.ATSTargetIterations ;
        else
            TimeStepDownRatio=NaN;
        end


        fprintf(CtrlVar.fidlog,' Adaptive Time Stepping:  #Non-Lin Iterations over last %-i time steps: (max|mean|min)=(%-g|%-g|%-g). Target is %-i. \t TimeStepUpRatio=%-g \n ',...
            nItVector,max(ItVector),mean(ItVector),min(ItVector),CtrlVar.ATSTargetIterations,TimeStepUpRatio);

        if RunInfo.Forward.uvhIterations(CtrlVar.CurrentRunStepNumber-1)==666  % This is a special forced reduction whenever RunInfo.Forward.uvhIterations has been set to this value


            dtOut=dtIn/CtrlVar.ATStimeStepFactorDown;
            RunInfo.Forward.AdaptiveTimeSteppingResetCounter=0;
            fprintf(CtrlVar.fidlog,' ---------------- Adaptive Time Stepping: time step decreased from %-g to %-g \n ',dtIn,dtOut);


        elseif RunInfo.Forward.AdaptiveTimeSteppingResetCounter > 2 && RunInfo.Forward.uvhIterations(CtrlVar.CurrentRunStepNumber-1)>25

            % This is also a special case to cover the possibilty that there is a sudden
            % increase in the number of non-linear iterations, or if the initial time step
            % a the start of a run was set too large.
            dtOut=dtIn/CtrlVar.ATStimeStepFactorDown;
            dtOut=max(dtOut,CtrlVar.ATSdtMin) ;
            RunInfo.Forward.AdaptiveTimeSteppingResetCounter=0;
            if dtOut<dtIn
                fprintf(' ---------------- Adaptive Time Stepping: time step decreased from %-g to %-g \n ',dtIn,dtOut);
            end
        else

            % This is the more general case.

            if RunInfo.Forward.AdaptiveTimeSteppingResetCounter>CtrlVar.ATSintervalDown && ~isnan(TimeStepDownRatio)

                % Potentially decrease time step

                if all(ItVector(1:CtrlVar.ATSintervalDown) > (CtrlVar.ATSTargetIterations+2) )  ||  ( TimeStepDownRatio > 2 )
                    dtOut=dtIn/CtrlVar.ATStimeStepFactorDown;
                    dtOut=max(dtOut,CtrlVar.ATSdtMin) ;
                    RunInfo.Forward.AdaptiveTimeSteppingResetCounter=0;

                    if dtOut<dtIn
                        fprintf(' ---------------- Adaptive Time Stepping: time step decreased from %-g to %-g \n ',dtIn,dtOut)
                    end

                end
            end

            if RunInfo.Forward.AdaptiveTimeSteppingResetCounter>CtrlVar.ATSintervalUp  && ~isnan(TimeStepUpRatio)
                if  TimeStepUpRatio<1

                    % Potentially increase time step

                    dtOut=min(CtrlVar.ATSdtMax,dtIn*CtrlVar.ATStimeStepFactorUp);
                    RunInfo.Forward.AdaptiveTimeSteppingResetCounter=0;



               


                end
            end
        end
    end

end

 dtOut=round(dtOut,2,"significant") ;



if  CtrlVar.UaRunType=="-uv-h-"  || CtrlVar.EnforceCFL    % If in semi-implicit step, make sure not to violate CFL condition

    dtcritical=CalcCFLdt2D(UserVar,RunInfo,CtrlVar,MUA,F) ;

    dtcritical=round(dtcritical,2,"significant") ;
    if ~isnan(dtcritical)
    
        nFactorSafety=2;

        if dtOut>dtcritical/nFactorSafety

            dtOut=dtcritical/nFactorSafety ;
            dtOut=round(dtOut,2,"significant") ;
           % fprintf('AdaptiveTimeStepping: dt > dt (CFL) and therefore dt reduced to %f \n',dtOut)

        end


        dtOut=min(dtcritical/nFactorSafety,dtOut*1.2) ;  % don't increase time step by more than 20%
        % fprintf('AdaptiveTimeStepping: dt=dtCFL/%i=%f \n',nFactorSafety,dtOut)

    end
end


RunInfo.Forward.dtRestart=dtOut ;  % Create a copy of dtOut before final modifications related to plot times and end times.
% This is the dt to be used in further restart runs

%% dtOut has now been set, but I need to see if the user wants outputs/plots at given time intervals and
% if I am possibly overstepping one of those intervals.
%
dtOutCopy=dtOut;  % keep a copy of dtOut to be able to revert to previous time step
% after this adjustment




%
if CtrlVar.DefineOutputsDt>0

    % Unlikely to be the case, but make sure dtOut is equal or smaller than Dt
    if dtOut>CtrlVar.DefineOutputsDt
        dtOut=CtrlVar.DefineOutputsDt;
    end

    %temp=dtOut;
    %dtOut=NoOverStepping(CtrlVar,time,dtOutCopy,CtrlVar.DefineOutputsDt);
    dtNoOverStepping=NoOverStepping(CtrlVar,time,dtOut,CtrlVar.DefineOutputsDt);
    % using dtNoOverStepping ensures that next output interval is not stepped over
    % dtNoOverStepping might be smaller that dtOut to avoid overstepping, but never smaller
    %
    %
    %   dtNoOverStepping is the dt I need to not overstep the Dt user requirements
    %  But do I need to change dt?  What if dtNoOverStepping is 'small'
    %
    % I use:
    %
    %   (ReminderFraction(CtrlVar.time,CtrlVar.DefineOutputsDt) < (CtrlVar.dt/(10*CtrlVar.DefineOutputsDt)))
    %
    % to determine if I need an output file.  So if an output file is created for both, or neither, dtOut and dtNoOverStepping,
    % there is no need to change dt
    %

    if dtNoOverStepping<dtOut  % so I may need to reduce dtOut, but maybe dtOut is only be a bit smaller than dtNoOverStepping


        % would I call DefineOutputs anyhow, even if I don't change dt? If so, then there is no need to change dtOut
        % and dtOut and dtNoOverStepping are too similar for the difference to matter.
        T1=ReminderFraction(CtrlVar.time+dtOut,CtrlVar.DefineOutputsDt) < (dtOut/(100*CtrlVar.DefineOutputsDt)) ;
        T2=ReminderFraction(CtrlVar.time+dtNoOverStepping,CtrlVar.DefineOutputsDt) < (dtNoOverStepping/(100*CtrlVar.DefineOutputsDt)) ;

        if T2  && ~T1
            dtOut=dtNoOverStepping;
            fprintf(CtrlVar.fidlog,' Adaptive Time Stepping: dt modified to accommodate user output requirements and set to %-g \n ',dtOut);
            RunInfo.Forward.AdaptiveTimeSteppingTimeStepModifiedForOutputs=1;
        elseif ~T2
            fprinft("I think this should not happen")
        else
            fprintf(" [dtNoOverStepping dtOut]=[%f %f]\n",dtNoOverStepping,dtOut);
        end
    elseif dtNoOverStepping>dtOut
        % it is possible that dtNoOverStepping is larger than dtOut. This happens if the remaining time to next output interval is
        % small, in which case dtNoOverStepping is increased for the upcoming time step to reach that output time.
        dtOut=dtNoOverStepping;
    end
end


if CtrlVar.ATSdtMax <= dtOut
    dtOut=CtrlVar.ATSdtMax ;
    fprintf(CtrlVar.fidlog,' ---------------- Adaptive Time Stepping: time step has reaced max allowed automated time step of %-g and is therefore not increased further \n ',CtrlVar.ATSdtMax);
else
    if dtOut>dtIn
        fprintf(CtrlVar.fidlog,' ---------------- Adaptive Time Stepping: time step increased from %-g to %-g \n ',dtIn,dtOut);
    end
end
%% make sure that run time does not exceed total run time as defined by user
% also check if current time is very close to total time, in which case there
% is no need to change the time step
if time+dtOut>CtrlVar.TotalTime && abs(time-CtrlVar.TotalTime)>100*eps

    dtOutOld=dtOut;
    dtOut=CtrlVar.TotalTime-time;

    if dtOutOld ~= dtOut
        fprintf(CtrlVar.fidlog,' Adaptive Time Stepping: dt modified to %-g to give a correct total run time of %-g \n ',dtOut,CtrlVar.TotalTime);
    end
end


%%

if dtOutCopy~=dtOut
    dtNotUserAdjusted=dtOutCopy;
else
    dtNotUserAdjusted=[];
end

dtOutLast=dtOut;


dtRatio=dtOut/dtIn;

%%
if dtOut==0
    save TestSave
    error('dtOut is zero')
end




end



