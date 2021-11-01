function [stop,Outs] = fminuncOutfun(p,optimValues,state)

persistent pOuts iCounter X StoreSolution

if nargin>0
    if isfield(p,'Inverse')
        if isfield(p.Inverse,'StoreSolutionAtEachIteration')
            if p.Inverse.StoreSolutionAtEachIteration
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
    pOuts.GradNorm=zeros(1000,1)+NaN;
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
        pOuts.GradNorm(iCounter)=norm(optimValues.gradient)/sqrt(numel(optimValues.gradient));
    end
    if StoreSolution
        X{iCounter}=p;
    else
        X=[];
    end
    iCounter=iCounter+1;
    
else
    Outs=pOuts;
    Outs.iteration=Outs.iteration(~isnan(Outs.iteration));
    Outs.fval=Outs.fval(~isnan(Outs.fval));
    Outs.StepSize=Outs.StepSize(~isnan(Outs.StepSize));
    Outs.GradNorm=Outs.GradNorm(~isnan(Outs.GradNorm));
    Outs.p=X;
end




end