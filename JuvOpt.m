function [F,J,RunInfo]=JuvOpt(UserVar,RunInfo,CtrlVar,MUA,BCs,F)


Ex="Quad";   % Ex="";
Opt="Mat" ; Opt="Ua" ; 


A=[] ; b=[] ; uvlb=[] ; uvub=[] ;  nonlcon = [];
[Aeq,beq]=AssembleLuvSSTREAM(CtrlVar,MUA,BCs) ;


fext0=KRTFgeneralBCs(CtrlVar,MUA,F,true); % RHS with velocities set to zero, i.e. only external forces

CtrlVar.Ex=Ex;



if CtrlVar.Ex=="Quad"
    func=@(uv) JQuad(uv,UserVar,RunInfo,CtrlVar,MUA,F,fext0) ;
   

else
    func=@(uv) Juv(uv,UserVar,RunInfo,CtrlVar,MUA,F,fext0,Aeq,beq) ;
    
end
 
uv=[] ; lambda=[];
Hess=Huv(uv,lambda,UserVar,RunInfo,CtrlVar,MUA,F) ;


neval=10;
nit=10;

% fmincon

optionsFminconInterior = optimoptions('fmincon',...
    'Algorithm','interior-point',...
    'CheckGradients',false,...
    'ConstraintTolerance',1e-6,...
    'HonorBounds',true,...
    'Diagnostics','on',...
    'DiffMaxChange',Inf,...
    'DiffMinChange',0,...
    'Display','iter-detailed',...
    'FunValCheck','off',...
    'MaxFunctionEvaluations',neval,...
    'MaxIterations',nit,...,...
    'OptimalityTolerance',1e-20,...
    'OutputFcn',@fminuncOutfun,...
    'PlotFcn',{@optimplotlogfval,@optimplotstepsize},...
    'HessianFcn',@Huv,...
    'SpecifyConstraintGradient',false,...
    'SpecifyObjectiveGradient',true,...
    'SubproblemAlgorithm','cg');  % here the options are 'gc'


optionsFminconTrust = optimoptions('fmincon',...
    'Algorithm','trust-region-reflective',...
    'CheckGradients',false,...
    'ConstraintTolerance',1e-6,...
    'HonorBounds',true,...
    'Diagnostics','on',...
    'DiffMaxChange',Inf,...
    'DiffMinChange',0,...
    'Display','iter-detailed',...
    'FunValCheck','off',...
    'MaxFunctionEvaluations',neval,...
    'MaxIterations',nit,...,...
    'OptimalityTolerance',1e-20,...
    'OutputFcn',@fminuncOutfun,...
    'PlotFcn',{@optimplotlogfval,@optimplotstepsize},...
    'StepTolerance',1e-30,...
    'FunctionTolerance',1e-30,...
    'UseParallel',true,...
    'HessianFcn','objective',...
    'HessianMultiplyFcn',[],...
    'SpecifyConstraintGradient',false,...
    'SpecifyObjectiveGradient',true,...
    'SubproblemAlgorithm','cg');  % here the options are 'gc'

optionsFminunc = optimoptions('fminunc',...
    'Algorithm','trust-region',...
    'CheckGradients',false,...
    'Diagnostics','on',...
    'DiffMaxChange',Inf,...
    'DiffMinChange',0,...
    'Display','iter-detailed',...
    'FunValCheck','off',...
    'MaxFunctionEvaluations',neval,...
    'MaxIterations',nit,...,...
    'OptimalityTolerance',1e-20,...
    'OutputFcn',@fminuncOutfun,...
    'PlotFcn',{@optimplotlogfval,@optimplotstepsize},...
    'StepTolerance',1e-30,...
    'FunctionTolerance',1e-30,...
    'UseParallel',true,...
    'HessianFcn','objective',...
    'HessianMultiplyFcn',[],...
    'SpecifyObjectiveGradient',true,...
    'SubproblemAlgorithm','factorization');  % here the options are 'gc'



if Ex=="Quad"
    Aeq=[]; beq=[] ;
end

uv0=[F.ub ; F.vb];

if Opt=="Mat"

    % [uv,J,exitflag,output] = fmincon(func,uv0,A,b,Aeq,beq,uvlb,uvub,nonlcon,optionsFminconInterior);
    [uv,J,exitflag,output] = fmincon(func,uv0,A,b,Aeq,beq,uvlb,uvub,nonlcon,optionsFminconTrust);
    % [uv,J,exitflag,output] = fminunc(func,uv0,optionsFminunc);

else

    [uv,J,exitflag,output] = fminconUa(func,uv0,A,b,Aeq,beq,uvlb,uvub,nonlcon,CtrlVar);

end


N=MUA.Nnodes;

F.ub=uv(1:N);
F.vb=uv(N+1:2*N);


end



