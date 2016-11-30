function [RunInfo,dtOut,dtRatio]=AdaptiveTimeStepping(RunInfo,CtrlVar,time,dtIn)
            
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
    %  -time step is not increased further than the target time step CtrlVar.ATStimeStepTarget
    %  -time step is adjusted so that total simulation time does not exceed CtrlVar.TotalTime
    % 
    %
    %
    
    persistent ItVector icount dtNotUserAdjusted dtOutLast
    
    % potentially dt was previously adjusted for plotting/saving purposes
    % if so then dtNotUserAdjusted is the previous unmodified time step
    % and I do want to revert to that time step. I do so only
    % however if dt has not been changed outside of function
    if ~isempty(dtOutLast) 
        if dtIn==dtOutLast  % dt not modified outside of function
            if ~isempty(dtNotUserAdjusted)   
                dtIn=dtNotUserAdjusted ; 
            end
        end
    end
    
    dtOut=dtIn ;
    
    
    %%
    if CtrlVar.AdaptiveTimeStepping && ~isempty(RunInfo) 
        
        if isempty(ItVector) ; ItVector=zeros(max(CtrlVar.ATSintervalDown,CtrlVar.ATSintervalUp),1)+1e10; end
        if isempty(icount) ; icount=0 ; end
        
        
        
        icount=icount+1;
        ItVector(2:end)=ItVector(1:end-1) ; ItVector(1)=RunInfo.Iterations;
        nItVector=numel(find(ItVector<1e10)); 
        TimeStepUpRatio=max(ItVector(1:CtrlVar.ATSintervalUp))/CtrlVar.ATSTargetIterations ;
        
        fprintf(CtrlVar.fidlog,' Adaptive Time Stepping:  #Non-Lin Iterations over last %-i time steps: (max|mean|min)=(%-g|%-g|%-g). Target is %-i. \t TimeStepUpRatio=%-g \n ',...
            nItVector,max(ItVector),mean(ItVector),min(ItVector),CtrlVar.ATSTargetIterations,TimeStepUpRatio);
        
        if icount>2 && RunInfo.Iterations>25
            icount=0;
            dtOut=dtIn/CtrlVar.ATStimeStepFactorDown;
            fprintf(CtrlVar.fidlog,' ---------------- Adaptive Time Stepping: time step decreased from %-g to %-g \n ',dtIn,dtOut);
        else
            
            if icount>CtrlVar.ATSintervalDown && nItVector >= CtrlVar.ATSintervalDown
                %if mean(ItVector(1:CtrlVar.ATSintervalDown)) > 5
                if all(ItVector(1:CtrlVar.ATSintervalDown) > 10 )
                    dtOut=dtIn/CtrlVar.ATStimeStepFactorDown; icount=0;
                    
                    
                    fprintf(CtrlVar.fidlog,' ---------------- Adaptive Time Stepping: time step decreased from %-g to %-g \n ',dtIn,dtOut);
                end
            end
            
            if icount>CtrlVar.ATSintervalUp && nItVector >= CtrlVar.ATSintervalUp
                if  TimeStepUpRatio<1
                    dtOut=min(CtrlVar.ATStimeStepTarget,dtIn*CtrlVar.ATStimeStepFactorUp); icount=0;

                    
                    if  CtrlVar.UaOutputsDt>0  
                        
                        Fraction=dtOut/CtrlVar.UaOutputsDt;
                        if Fraction>=1  % Make sure dt is not larger than the interval between UaOutputs
                            dtOut=CtrlVar.UaOutputsDt;
                        elseif Fraction>0.1   % if dt is greater than 10% of UaOutputs interval, round dt
                                              % so that it is an interger multiple of UaOutputsDt
                            fprintf('Adaptive Time Stepping dtout=%f \n ',dtOut);
                            dtOut=CtrlVar.UaOutputsDt/RoundNumber(CtrlVar.UaOutputsDt/dtOut,1);
                            fprintf('Adaptive Time Stepping dtout=%f \n ',dtOut);
                        end
                    end
                    
                    
                    
                    if CtrlVar.ATStimeStepTarget <= dtOut
                        dtOut=CtrlVar.ATStimeStepTarget ;
                        fprintf(CtrlVar.fidlog,' ---------------- Adaptive Time Stepping: time step has reached target time step of %-g and is therefore not increased further \n ',CtrlVar.ATStimeStepTarget);
                    else
                        fprintf(CtrlVar.fidlog,' ---------------- Adaptive Time Stepping: time step increased from %-g to %-g \n ',dtIn,dtOut);
                    end
                    
                    
                end
            end
        end
    end
    
    %% dtOut has now been set, but I need to see if the user wants outputs/plots at given time intervals and
    % if I am possibly overstepping one of those intervals.
    %
    dtOutCopy=dtOut;  % keep a copy of dtOut to be able to revert to previous time step
    % after this adjustment
    
    
    if CtrlVar.WriteDumpFile && CtrlVar.WriteDumpFileTimeInterval>0
        temp=dtOut;
        dtOut=NoOverStepping(CtrlVar,time,dtOut,CtrlVar.WriteDumpFileTimeInterval);
        if abs(temp-dtOut)>100*eps
            fprintf(CtrlVar.fidlog,' Adaptive Time Stepping: dt modified to accomondate output file requirements and set to %-g \n ',dtOut);
        end
    end
    
    %
    if CtrlVar.UaOutputsDt>0
        temp=dtOut;
        dtOut=NoOverStepping(CtrlVar,time,dtOutCopy,CtrlVar.UaOutputsDt);
        if abs(temp-dtOut)>100*eps
            fprintf(CtrlVar.fidlog,' Adaptive Time Stepping: dt modified to accomondate user output requirements and set to %-g \n ',dtOut);
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



