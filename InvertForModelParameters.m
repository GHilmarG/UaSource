
function [UserVar,F,l,InvFinalValues,RunInfo]=...
    InvertForModelParameters(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

narginchk(11,11)

InvFinalValues=InvStartValues;

if isempty(CtrlVar.Inverse.InitialLineSearchStepSize) ||  CtrlVar.Inverse.InitialLineSearchStepSize==0
    CtrlVar.Inverse.InitialLineSearchStepSize=InvStartValues.SearchStepSize;
end


%% Consider adding some test of input variables here, maybe see if all gs,ga fields are defined

%% Define inverse parameters and anonymous function returning objective function, directional derivative, and Hessian
%


F=InvStartValues2F(CtrlVar,MUA,F,InvStartValues,Priors,Meas) ;
[F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);


[F.GF,GLgeo,GLnodes,GLele]=IceSheetIceShelves(CtrlVar,MUA,F.GF) ;

% p is the vector of the control variables, currenty p=[A,b,C]
% with A, b or C here only being nonempty when inverted for, 
% This mapping between A, b and C into the control variable is done by F2p

% Make sure initial point is feasible
F.AGlen=kk_proj(F.AGlen,F.AGlenmax,F.AGlenmin) ;
F.C=kk_proj(F.C,F.Cmax,F.Cmin) ;

[p0,plb,pub]=F2p(CtrlVar,MUA,F); 


CtrlVar.Inverse.ResetPersistentVariables=1;
[J0,dJdp,Hessian,JGHouts,F,RunInfo]=JGH(p0,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
CtrlVar.Inverse.ResetPersistentVariables=0;
% The parameters passed in the anonymous function are those that exist at the time the anonymous function is created.


CtrlVar.WriteRunInfoFile=0;
func=@(p) JGH(p,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);

Hfunc=@(p,lambda) HessianAC(p,lambda,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);

fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \n \n',J0,JGHouts.MisfitOuts.I,JGHouts.RegOuts.R,norm(dJdp))

dJdpTest=[];

%%

if CtrlVar.Inverse.TestAdjoint.isTrue
    
   
    % Get the gradient using the adjoint method
    [J,dJdp,Hessian,JGHouts]=func(p0);
    
    
    NA=numel(InvStartValues.AGlen);  % Number of A parameters to invert for
    
    if isempty(CtrlVar.Inverse.TestAdjoint.iRange)
        iRange=1:NA;
    else
        iRange=CtrlVar.Inverse.TestAdjoint.iRange;
    end
    
    switch strlength(CtrlVar.Inverse.InvertForField)
        
        case 2
            iRange=[iRange(:);iRange(:)+NA];
        case 3
            iRange=[iRange(:);iRange(:)+NA;iRange(:)+2*NA];
    end
    
    I=(iRange>=1) & (iRange <= numel(p0));
    iRange=iRange(I);
    
    % calc brute force gradient

    
    dJdpTest = CalcBruteForceGradient(func,p0,CtrlVar,iRange);

    filename=CtrlVar.Experiment+"BruteForceGradient";
    fprintf('BruteForceGradient save in the file : %s \n',filename)
    save(filename,'CtrlVar','UserVar','MUA','F','dJdpTest','iRange')
    
    
else
    
    
    %%
    
    if contains(CtrlVar.Inverse.MinimisationMethod,"Ua")
        
        [p,UserVar,RunInfo]=UaOptimisation(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub);
        
        
    elseif contains(CtrlVar.Inverse.MinimisationMethod,"Matlab")
        
        clear fminconOutputFunction fminconHessianFcn fminuncOutfun
        
        [p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub,Hfunc);
        
        
    else
        
        fprintf(' CtrlVar.Inverse.MinimisationMethod has the value %s \n',CtrlVar.Inverse.MinimisationMethod)
        fprintf(' but can only have the values ''MatlabOptimization'' or ''UaOptimization''\n')
        error('what case? ')
    end
    
    F=p2F(CtrlVar,MUA,p,F,Meas,Priors);
    [J,dJdp,Hessian,JGHouts,F]=JGH(p,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
    fprintf('\n +++++++++++ At end of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \n \n',J,JGHouts.MisfitOuts.I,JGHouts.RegOuts.R,norm(dJdp))
    
    
end

% Put RAa, RAs, RCa, RCs in InvFinalValues
InvFinalValues=Vars2InvValues(CtrlVar,F,InvFinalValues,J,dJdp,JGHouts,RunInfo,dJdpTest); 


end
