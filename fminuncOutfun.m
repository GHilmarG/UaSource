function [stop,Outs] = fminuncOutfun(x,optimValues,state)

persistent pOuts iCounter X StoreSolution

if nargin>0
    if isfield(x,'Inverse')
        if isfield(x.Inverse,'StoreSolutionAtEachIteration')
            if x.Inverse.StoreSolutionAtEachIteration
                StoreSolution=1;
                return
            end
        end
    end
end


if isempty(pOuts)
    pOuts.fval=zeros(1000,1)+NaN;
    pOuts.iteration=zeros(1000,1)+NaN;
    pOuts.StepSize=zeros(1000,1)+NaN;
    iCounter=1;
end

stop=false;

if nargin>0
    
    pOuts.fval(iCounter)=optimValues.fval;
    pOuts.iteration(iCounter)=optimValues.iteration;
    if isempty(optimValues.stepsize)
        pOuts.StepSize(iCounter)=NaN;
    else
        pOuts.StepSize(iCounter)=optimValues.stepsize;
    end
    if StoreSolution
        X{iCounter}=x;
    else
        X=[];
    end
    iCounter=iCounter+1;
    
else
    Outs=pOuts;
    Outs.iteration=Outs.iteration(~isnan(Outs.iteration));
    Outs.fval=Outs.fval(~isnan(Outs.fval));
    Outs.StepSize=Outs.StepSize(~isnan(Outs.StepSize));
    Outs.p=X;
end





end