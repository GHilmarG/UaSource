function [stop,info] = fminconOutputFunction(x,optimValues,state)


persistent Info


if ~isempty(optimValues)
    if isempty(Info)
        Info.fval=[] ; Info.stepsize=[]; Info.iteration=[];
    end
    
    
    if optimValues.iteration>0
        Info.fval(optimValues.iteration)=optimValues.fval;
        Info.stepsize(optimValues.iteration)= optimValues.stepsize;
        
    end
    
end

stop=0;

info=Info;


end