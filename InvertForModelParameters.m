
function [UserVar,F,l,InvFinalValues,RunInfo]=...
    InvertForModelParameters(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

InvFinalValues=InvStartValues;

if isempty(CtrlVar.Inverse.InitialLineSearchStepSize) ||  CtrlVar.Inverse.InitialLineSearchStepSize==0
    CtrlVar.Inverse.InitialLineSearchStepSize=InvStartValues.SearchStepSize;
end


%% Did user define regularization parameters in `DefineInputForInverseRun'?
% If so, then use those values

if ~isempty(Priors.Regularize)
    
    if ~isempty(Priors.Regularize.logC.gs)
        CtrlVar.Inverse.Regularize.logC.gs=Priors.Regularize.logC.gs;
    end
    
    if ~isempty(Priors.Regularize.logC.ga)
        CtrlVar.Inverse.Regularize.logC.ga=Priors.Regularize.logC.ga;
    end
    
    
    if ~isempty(Priors.Regularize.C.gs)
        CtrlVar.Inverse.Regularize.C.gs=Priors.Regularize.C.gs;
    end
    
    if ~isempty(Priors.Regularize.C.ga)
        CtrlVar.Inverse.Regularize.C.ga=Priors.Regularize.C.ga;
    end
    
    
    if ~isempty(Priors.Regularize.logAGlen.gs)
        CtrlVar.Inverse.Regularize.logAGlen.gs=Priors.Regularize.logAGlen.gs;
    end
    
    if ~isempty(Priors.Regularize.logAGlen.ga)
        CtrlVar.Inverse.Regularize.logAGlen.ga=Priors.Regularize.logAGlen.ga;
    end
    
    
    if ~isempty(Priors.Regularize.AGlen.gs)
        CtrlVar.Inverse.Regularize.AGlen.gs=Priors.Regularize.AGlen.gs;
    end
    
    if ~isempty(Priors.Regularize.AGlen.ga)
        CtrlVar.Inverse.Regularize.AGlen.ga=Priors.Regularize.AGlen.ga;
    end
    
end

%% Define inverse parameters and anonymous function returning objective function, directional derivative, and Hessian
%




[p0,plb,pub]=InvValues2p(CtrlVar,InvStartValues); 



CtrlVar.Inverse.ResetPersistentVariables=1;
[J0,dJdp,Hessian,JGHouts,F]=JGH(p0,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);

% The parameters passed in the anonymous function are those that exist at the time the anonymous function is created.



func=@(p) JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);

fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \n \n',J0,JGHouts.MisfitOuts.I,JGHouts.RegOuts.R,norm(dJdp))


CtrlVar.Inverse.ResetPersistentVariables=0;
%%

if CtrlVar.Inverse.TestAdjoint.isTrue
    
    
    % Get the gradient using the adjoint method
    [J,dJdp,Hessian]=func(p0);
    
    
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
            
            [p,RunInfo]=UaOptimisation(CtrlVar,func,p0,plb,pub,RunInfo);
            
            
            
        case 'MatlabOptimization'
            
            clear fminconOutputFunction fminconHessianFcn fminuncOutfun
            
            [p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(CtrlVar,func,p0,plb,pub,RunInfo);
            
            
            
        otherwise
            
            fprintf(' CtrlVar.Inverse.MinimisationMethod has the value %s \n',CtrlVar.Inverse.MinimisationMethod)
            fprintf(' but can only have the values ''MatlabOptimization'' or ''UaOptimization''\n')
            error('what case? ')
    end
    
    
    
    %     [InvFinalValues.C,iU,iL]=kk_proj(InvFinalValues.C,CtrlVar.Cmax,CtrlVar.Cmin);
    %     [InvFinalValues.AGlen,iU,iL]=kk_proj(InvFinalValues.AGlen,CtrlVar.AGlenmax,CtrlVar.AGlenmin);
    %
    

    [J,dJdp,Hessian,JGHouts,F]=JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
    fprintf('\n +++++++++++ At end of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \n \n',J,JGHouts.MisfitOuts.I,JGHouts.RegOuts.R,norm(dJdp))
    
    
    NA=numel(InvStartValues.AGlen);
    NC=numel(InvStartValues.C);
    InvFinalValues=p2InvValues(CtrlVar,p,InvFinalValues,NA,NC);

  
    
end

%% Values
InvFinalValues.J=J;
InvFinalValues.I=JGHouts.MisfitOuts.I;
InvFinalValues.R=JGHouts.RegOuts.R;
InvFinalValues.RAGlen=JGHouts.RegOuts.RAGlen;
InvFinalValues.RC=JGHouts.RegOuts.RC;
if isprop(InvFinalValues,'uAdjoint')
    InvFinalValues.uAdjoint=JGHouts.MisfitOuts.uAdjoint;
    InvFinalValues.vAdjoint=JGHouts.MisfitOuts.vAdjoint;
end


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
