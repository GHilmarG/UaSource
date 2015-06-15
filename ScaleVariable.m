function out=ScaleVariable(from,to,in)
    
    % scales input variable `in' so that on output its range is between [from,to]
    %
    
    if numel(in)==1
        fprintf('can not scale variable if only one value given \n')
        out=NaN;
        return
        
    elseif max(in)==min(in)

        out=NaN;
        
    else
        
        out=from+(to-from).*(in-min(in))./(max(in)-min(in));
        
    end
end
