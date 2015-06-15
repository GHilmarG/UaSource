function [InvStartValues,Priors,Meas,BCsAdjoint]=GetInputsForInverseRun(Experiment,CtrlVar,MUA,BCs,time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF)


listing=dir('DefineInputsForInverseRun.m') ;

BCsAdjoint=BoundaryConditions;
Meas=Measurements;
Priors=PriorProbabilityDistribution;
InvStartValues=InversionStartValues;

if numel(listing)==0
    
    
    [s_prior,b_prior,S_prior,B_prior,Priors.AGlen,Priors.C,Priors.n,Priors.m,rho_prior,rhow_prior,...
        Priors.CovAGlen,Priors.CovC,...
        lxFixedNode,lyFixedNode,lxFixedValue,lyFixedValue,lxtiedA,lxtiedB,lytiedA,lytiedB,...
        s_meas,Meas.us,Meas.vs,Meas.ws,b_meas,B_meas,...
        usError,vsError,wsError]=...
        DefineInverseModellingVariables(Experiment,CtrlVar,MUA,time,...
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
       
else
    
    [InvStartValues,Priors,Meas,BCsAdjoint]=DefineInputsForInverseRun(Experiment,CtrlVar,MUA,BCs,InvStartValues,Priors,Meas,BCsAdjoint,time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF);
    
    
end


if isempty(InvStartValues.AGlen) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.AGlen is empty') ; end
if isempty(InvStartValues.C) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.C is empty') ; end
if isempty(InvStartValues.n) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.n is empty') ; end
if isempty(InvStartValues.m) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.m is empty') ; end


end

