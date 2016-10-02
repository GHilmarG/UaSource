function [UserVar,InvStartValues,Priors,Meas,BCsAdjoint]=GetInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF,g,alpha,ub,vb,ud,vd,l)



BCsAdjoint=BoundaryConditions;
Meas=Measurements;
Priors=PriorProbabilityDistribution;
InvStartValues=InversionStartValues;

if exist(fullfile(cd,'DefineInputsForInverseRun.m'),'file')
    
    
    N=nargout('DefineInputsForInverseRun');
    
    switch N
        
        case 4
            
            [InvStartValues,Priors,Meas,BCsAdjoint]=DefineInputsForInverseRun(CtrlVar.UserVar,CtrlVar,MUA,BCs,InvStartValues,Priors,Meas,BCsAdjoint,time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF);
            
        case 5
            
            [UserVar,InvStartValues,Priors,Meas,BCsAdjoint]=DefineInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,InvStartValues,Priors,Meas,BCsAdjoint,time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF,g,alpha,ub,vb,ud,vd,l);
    
    end
    
    BCsAdjoint=CreatePlausibleBCsForAdjointProblem(BCs,BCsAdjoint);
    
    
else
    
    fprintf('\n-----------------------------------------------------------------------')
    fprintf(' Found ''DefineInverseModellingVariables.m'' in local run directory.\n')
    fprintf(' That m-file will be used to define inputs for inversion.\n')
    fprintf(' However, this is no longer the recomended approach.\n')
    fprintf(' Suggest switching to using ''DefineInputsForInverseRun.m'' \n')
    fprintf('-----------------------------------------------------------------------\n')
    
    [s_prior,b_prior,S_prior,B_prior,Priors.AGlen,Priors.C,Priors.n,Priors.m,Priors.rho,Priors.rhow,...
        Priors.CovAGlen,Priors.CovC,...
        lxFixedNode,lyFixedNode,lxFixedValue,lyFixedValue,lxtiedA,lxtiedB,lytiedA,lytiedB,...
        s_meas,Meas.us,Meas.vs,Meas.ws,b_meas,B_meas,...
        usError,vsError,wsError]=...
        DefineInverseModellingVariables(UserVar,CtrlVar,MUA,time,...
        BCs.ubFixedNode,BCs.ubFixedValue,BCs.vbFixedNode,BCs.vbFixedValue,BCs.ubTiedNodeA,BCs.ubTiedNodeB,BCs.vbTiedNodeA,BCs.vbTiedNodeB,...
        AGlen,C,n,m,s,b,S,B,rho,rhow,GF);
    
    BCsAdjoint.ubFixedNode=lxFixedNode;
    BCsAdjoint.vbFixedNode=lyFixedNode;
    BCsAdjoint.ubFixedValue=lxFixedValue;
    BCsAdjoint.vbFixedValue=lyFixedValue;
    BCsAdjoint.ubTiedNodeA=lxtiedA;
    BCsAdjoint.ubTiedNodeB=lxtiedB;
    BCsAdjoint.vbTiedNodeA=lytiedA;
    BCsAdjoint.vbTiedNodeB=lytiedB;
    
    
    %InvStartValues.B=Priors.B;
    InvStartValues.m=Priors.m;
    InvStartValues.n=Priors.n;
    InvStartValues.C=Priors.C;
    InvStartValues.AGlen=Priors.AGlen;
    
    
    Meas.usCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,usError.^2,MUA.Nnodes,MUA.Nnodes);
    Meas.vsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,vsError.^2,MUA.Nnodes,MUA.Nnodes);
    Meas.wsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,wsError.^2,MUA.Nnodes,MUA.Nnodes);
    
    
    
end

%% Test inputs
Meas=TestMeas(CtrlVar,MUA,Meas);

[InvStartValues.AGlen,InvStartValues.n]=TestAGlenInputValues(CtrlVar,MUA,InvStartValues.AGlen,InvStartValues.n);
[Priors.AGlen,Priors.n]=TestAGlenInputValues(CtrlVar,MUA,Priors.AGlen,Priors.n);

[InvStartValues.C,InvStartValues.m]=TestSlipperinessInputValues(CtrlVar,MUA,InvStartValues.C,InvStartValues.m);
[Priors.C,Priors.m]=TestSlipperinessInputValues(CtrlVar,MUA,Priors.C,Priors.m);

[Priors.rho,Priors.rhow]=TestDensityInputValues(CtrlVar,MUA,Priors.rho,Priors.rhow);

isCorrectDimentions=DoPriorsHaveCorrectDimentions(CtrlVar,MUA,Priors);
if ~ isCorrectDimentions
    fprintf(' Priors do not have right dimentions at restart. \n')
    fprintf(' Modify DefineInputsForInverseRun to ensure that dimentions are correct.\n')
    error('Ua:GetInputForInverseRun:incorrectdimentisons','incorrect dimentions')
    
end

if isempty(InvStartValues.AGlen) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.AGlen is empty') ; end
if isempty(InvStartValues.C) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.C is empty') ; end
if isempty(InvStartValues.n) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.n is empty') ; end
if isempty(InvStartValues.m) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.m is empty') ; end


end

























