function [UserVar,Cest,AGlenEst,Info,ub,vb,ud,vd,l,xAdjoint,yAdjoint,gammaAdjoint]=InversionUsingMatlabOptimizationToolbox(...
    UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info)

%%

CtrlVar.OnlyGetPersistenValues=0;
xAdjoint=[] ;yAdjoint=[]; gammaAdjoint=[];


AGlen=InvStartValues.AGlen;
C=InvStartValues.C;

AGlenEst=AGlen;
Cest=C;


n=InvStartValues.n;
m=InvStartValues.m;

switch CtrlVar.AdjointGrad
    
    case 'A'
        
        func=@(AGlen) CostFunctionValueAndGradient(...
            UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF,BCsAdjoint,Priors,Meas);
        
        lowerb=zeros(length(AGlen),1)+CtrlVar.AGlenmin-eps;
        upperb=zeros(length(AGlen),1)+CtrlVar.AGlenmax+eps;
        TypicalX=mean(AGlen)+zeros(numel(AGlen),1); 
        TolCon=CtrlVar.AGlenmin/10;
        InitTrustRegionRadius=median(C)/100;
        x0=InvStartValues.AGlen;
        
    case 'C'
        
        
        func=@(C) CostFunctionValueAndGradient(...
            UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF,BCsAdjoint,Priors,Meas);
        
        lowerb=zeros(length(C),1)+CtrlVar.Cmin*2;
        upperb=zeros(length(C),1)+CtrlVar.Cmax/2;
        TypicalX=mean(C)+zeros(numel(C),1); 
        TolCon=CtrlVar.Cmin/10;
        InitTrustRegionRadius=median(C)/100;
        x0=InvStartValues.C;
        
end


HessianFcn=@(x,lambda) fminconHessianFcn(x,lambda,MUA,AGlen,C,CtrlVar);

switch CtrlVar.AdjointMinimisationMethod
    
    
    
    case {'MatlabOptimizationToolbox:fmincon','MatlabOptimizationToolbox'}
        
        % fmincon
        options=optimoptions('fmincon',...
            'Algorithm','interior-point',...
            'CheckGradients',false,...
            'ConstraintTolerance',TolCon,...
            'Diagnostics','on',...
            'DiffMaxChange',Inf,...
            'DiffMinChange',0,...
            'Display','iter-detailed',...
            'FunValCheck','off',...
            'MaxFunctionEvaluations',100000,...
            'MaxIterations',CtrlVar.MaxAdjointIterations,...
            'OptimalityTolerance',1e-20,...
            'OutputFcn',@fminconOutputFunction,... %             'PlotFcn',@optimplotfval,...
            'StepTolerance',1e-20,...
            'FunctionTolerance',1e-06,...
            'TypicalX',TypicalX,...
            'UseParallel',true,...
            'HessianApproximation',{'lbfgs',30},'HessianFcn',HessianFcn,'HessianMultiplyFcn',[],...
            'ScaleProblem','none',...
            'InitTrustRegionRadius',InitTrustRegionRadius,...
            'SpecifyConstraintGradient',false,...
            'SpecifyObjectiveGradient',true,...
            'SubproblemAlgorithm','factorization');
        
        [x,Misfit1,exitflag,output,lambda,grad,hessian]=fmincon(func,x0,[],[],[],[],lowerb,upperb,[],options);
        
    case {'MatlabOptimizationToolbox:fminunc'}
        
        % fminunc
        options=optimoptions('fminunc',...
            'Algorithm','trust-region',...
            'CheckGradients',false,...
            'Display','iter-detailed',...
            'FiniteDifferenceStepSize',sqrt(eps),...
            'FiniteDifferenceType','forward',...
            'FunctionTolerance',1e-06,...
            'HessianFcn','objective',...
            'MaxFunctionEvaluations',100000,...
            'MaxIterations',CtrlVar.MaxAdjointIterations,...
            'OptimalityTolerance',1e-20,...
            'OutputFcn',@fminconOutputFunction,...
            'PlotFcn',@optimplotfval,...
            'SpecifyObjectiveGradient',true,...
            'StepTolerance',1e-20,...
            'SubproblemAlgorithm','cg',...
            'TypicalX',TypicalX)';
        
        [x,Misfit1,exitflag,output,grad,hessian]=fminunc(func,x0,options);
        
end


switch CtrlVar.AdjointGrad
    
    case 'A'
        AGlenEst=x;
    case 'C'
        Cest=x;
end



CtrlVar.OnlyGetPersistenValues=1;

[~,~,~,~,~,~,ub,vb,ud,vd,l,UserVar]=ObjFunctionValueAndGradient(...
    UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlenEst,Cest,n,m,alpha,rho,rhow,g,GF,BCsAdjoint,Priors,Meas);



[stop,fminconInfo] = fminconOutputFunction([],[],[]);

nI=size(Info.JoptVector,1);
nJ=numel(fminconInfo.fval);
if isempty(Info.JoptVector)
    Info.JoptVector=zeros(nJ,7)+NaN; 
    Info.JoptVector(:,1)=fminconInfo.fval;
else
    Info.JoptVector=[Info.JoptVector;zeros(nJ,7)+NaN];
    Info.JoptVector(nI+1:end,1)=fminconInfo.fval;
end




end
