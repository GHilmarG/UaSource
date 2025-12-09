function isCorrectDimensions=DoPriorsHaveCorrectDimensions(CtrlVar,MUA,Priors)

isCorrectDimensions=1;

% 
% if isempty(Priors.n)
%     
%     error(' Priors.n is empty.')
%     
% end
% 
% 
% if isempty(Priors.m)
%     
%     error(' Priors.m is empty.')
%     
% end

if isempty(Priors.C)
    
    error(' Priors.C is empty.')
    
end

if isempty(Priors.AGlen)
    
    error(' Priors.AGlen is empty.')
    
end

if strcmpi(CtrlVar.Inverse.Regularize.Field,'cov')
    
    if contains(upper(CtrlVar.Inverse.InvertFor),'C')
        
        [nCC,mCC]=size(Priors.CovC) ;
        
        if CtrlVar.CisElementBased
            isCorrectDimensions=(nCC==MUA.Nele && mCC==MUA.Nele);
        else
            isCorrectDimensions=(nCC==MUA.Nnodes && mCC==MUA.Nnodes);
        end
    end
    
    if contains(upper(CtrlVar.Inverse.InvertFor),'A')
        
        [nAA,mAA]=size(Priors.CovAGlen) ;
        
        if CtrlVar.AGlenisElementBased
            isCorrectDimensions=(nAA==MUA.Nele && mAA==MUA.Nele);
        else
            isCorrectDimensions=(nAA==MUA.Nnodes && mAA==MUA.Nnodes);
        end
        
    end
    
    
end
