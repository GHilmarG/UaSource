
function [UserVar,F,l,InvFinalValues,RunInfo]=...
    InvertForModelParameters(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

InvFinalValues=InvStartValues;

if isempty(CtrlVar.Inverse.InitialLineSearchStepSize) ||  CtrlVar.Inverse.InitialLineSearchStepSize==0
    CtrlVar.Inverse.InitialLineSearchStepSize=InvStartValues.SearchStepSize;
end

%% Define inverse parameters and anonymous function returning objective function, directional derivative, and Hessian
%


if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
    
    pA0=log10(InvStartValues.AGlen);
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
    
    pA0=InvStartValues.AGlen;
    
elseif ~contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
    
    pA0=[];
    
end


if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
    
    pC0=log10(InvStartValues.C);
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'c')
    
    pC0=InvStartValues.C;
    
elseif ~contains(lower(CtrlVar.Inverse.InvertFor),'c')
    
    pC0=[];
    
end

p0=[pA0;pC0];



%

CtrlVar.Inverse.ResetPersistentVariables=1;
[J0,dJdp,Hessian,JGHouts]=JGH(p0,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
CtrlVar.Inverse.ResetPersistentVariables=0;
% The parameters passed in the anonymous function are those that exist at the time the anonymous function is created.
func=@(p) JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);


%%

if CtrlVar.Inverse.TestAdjoint.isTrue
    
    
    % Get the gradient using the adjoint method
    [J,dJdp,H]=func(p0);
    
    
    % calc brute force gradient
    dJdpTest = CalcBruteForceGradient(func,p0,CtrlVar);
    InvFinalValues.dJdpTest=dJdpTest;
    
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'aglen') && contains(lower(CtrlVar.Inverse.InvertFor),'c')
        
        NA=numel(InvStartValues.AGlen);
        InvFinalValues.dJdAGlenTest=dJdpTest(1:NA);
        InvFinalValues.dJdCTest=dJdpTest(NA+1:end);
        
        
    elseif contains(lower(CtrlVar.Inverse.InvertFor),'aglen') && ~contains(lower(CtrlVar.Inverse.InvertFor),'c')
        
        InvFinalValues.dJdAGlenTest=dJdpTest;
        InvFinalValues.dJdCTest=[];
        
    elseif ~contains(lower(CtrlVar.Inverse.InvertFor),'aglen') && contains(lower(CtrlVar.Inverse.InvertFor),'c')
        
        InvFinalValues.dJdAGlenTest=[];
        InvFinalValues.dJdCTest=dJdpTest;
        
    else
        error('sfda')
    end
    
else
    
    
    %%
    
    switch CtrlVar.Inverse.MinimisationMethod
        
        case 'UaOptimization'
            
            [p,RunInfo]=UaOptimisation(CtrlVar,func,p0,RunInfo);
            
            
            
        case 'MatlabOptimization'
            
            clear fminconOutputFunction fminconHessianFcn fminuncOutfun
            
            [p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(CtrlVar,func,p0,RunInfo);
            
            
            
        otherwise
            
            fprintf(' CtrlVar.Inverse.MinimisationMethod has the value %s \n',CtrlVar.Inverse.MinimisationMethod)
            fprintf(' but can only have the values ''MatlabOptimization'' or ''UaOptimization''\n')
            error('what case? ')
    end
    
    
    
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'aglen') && contains(lower(CtrlVar.Inverse.InvertFor),'c')  % AC
        
        NA=numel(InvStartValues.AGlen);
        NC=numel(InvStartValues.C);
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            InvFinalValues.AGlen=10.^p(1:NA);
        else
            InvFinalValues.AGlen=p(1:NA);
        end
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            InvFinalValues.C=10.^p(NA+1:end);
        else
            InvFinalValues.C=p(NA+1:end);
        end
        
    elseif contains(lower(CtrlVar.Inverse.InvertFor),'aglen')   % A
        
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            InvFinalValues.AGlen=10.^p;
        else
            InvFinalValues.AGlen=p;
        end
        
        
    elseif contains(lower(CtrlVar.Inverse.InvertFor),'c')  % C
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
            InvFinalValues.C=10.^p;
        else
            InvFinalValues.C=p;
        end
    else
        fprintf(' CtrlVar.Inverse.InvertFor=%s \n',CtrlVar.Inverse.InvertFor)
        fprintf(' CtrlVar.Inverse.InvertFor does not have expected value.\n')
        error('InverForModelParameters:incorrect inputs')
    end
    
    
    
    [InvFinalValues.C,iU,iL]=kk_proj(InvFinalValues.C,CtrlVar.Cmax,CtrlVar.Cmin);
    [InvFinalValues.AGlen,iU,iL]=kk_proj(InvFinalValues.AGlen,CtrlVar.AGlenmax,CtrlVar.AGlenmin);
    
    
    [J,dJdp,Hessian,JGHouts,F]=JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
    
    
end

%% Values
InvFinalValues.J=J;
InvFinalValues.I=JGHouts.MisfitOuts.I;
InvFinalValues.R=JGHouts.RegOuts.R;
InvFinalValues.RAGlen=JGHouts.RegOuts.RAGlen;
InvFinalValues.RC=JGHouts.RegOuts.RC;

%% Gradients
InvFinalValues.dJdp=dJdp;

InvFinalValues.dIdp=JGHouts.dIdp;
InvFinalValues.dRdp=JGHouts.dRdp;

InvFinalValues.dJdAGlen=JGHouts.MisfitOuts.dIdAGlen+JGHouts.RegOuts.dRdAGlen;
InvFinalValues.dJdC=JGHouts.MisfitOuts.dIdC+JGHouts.RegOuts.dRdC;

%% These are of less interest, but can be added
%InvFinalValues.dIdAGlen=JGHouts.MisfitOuts.dIdAGlen;
%InvFinalValues.dIdC=JGHouts.MisfitOuts.dIdC;

%InvFinalValues.dRdAGlen=JGHouts.RegOuts.dRdAGlen;
%InvFinalValues.dRdC=JGHouts.RegOuts.dRdC;


InvFinalValues.SearchStepSize=RunInfo.Inverse.StepSize(end);



end
