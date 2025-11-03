

%%

% load("0220000-FR2020to2500-10km-uvh-Tri3-SlidWeertman-Duvh-MRIM6HadGEM2-abMask0A-P-BCVel-kH10000-TM0k2-Alim-Clim-Ca1-Cs100000-Aa1-As100000-VelITS120-BM3-SMB_RACHMO2k3_2km-.mat")


load("../UaTests/RadialIceCap/RestartRadialDome.mat","CtrlVarInRestartFile","MUA","F","RunInfo","l","BCs")

%load("../UaTests/MismipPlus/RestartMismipPlusRoughSteadyStateSolution.mat","CtrlVarInRestartFile","MUA","F","BCs","l","RunInfo")

load("../UaTests/GaussPeak/GaussPeakRestart.mat","CtrlVarInRestartFile","MUA","F","RunInfo","l","BCs")


CtrlVar=CtrlVarInRestartFile;
CtrlVar.FlowApproximation="SSTREAM";

%BCs.ubFixedNode=1; BCs.ubFixedValue=BCs.ubFixedNode*0;
%BCs.vbFixedNode=1; BCs.vbFixedValue=BCs.vbFixedNode*0;

UserVar=[];
RunInfo=UaRunInfo();
MUA=UpdateMUA(CtrlVar,MUA);
CtrlVar.Development.Pre2025uvAssembly=false;

CtrlVar.uvAssembly.ZeroFields=true;   CtrlVar.uvMatrixAssembly.Ronly=true ;
% fext0=KRTFgeneralBCs(CtrlVar,MUA,F);            % RHS with velocities set to zero, i.e. only external forces

[RunInfo,fext0]=uvMatrixAssembly(RunInfo,CtrlVar,MUA,F,BCs);   % RHS with velocities set to zero, i.e. only external forces
Normalisation=fext0'*fext0+1000*eps;

%F.ub=F.ub*0;
%F.vb=F.vb*0; 
tuv=tic;
[UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
tuv=toc(tuv);

Res=full(Ruv'*Ruv/Normalisation);
fprintf("uv Residual=%g \n",Res)



fun=@(x)  JuvLSQ(x,RunInfo,CtrlVar,MUA,F,BCs) ;

%F.ub(BCs.ubFixedNode)=BCs.ubFixedValue; F.vb(BCs.vbFixedNode)=BCs.vbFixedValue;
x0=[F.ub;F.vb];




options = optimoptions('lsqnonlin',...
    Display='iter-detailed',...
    SpecifyObjectiveGradient=true,...
    SpecifyConstraintGradient=true,...
    Algorithm='interior-point',...
    BarrierParamUpdate='predictor-corrector',...
    OptimalityTolerance=1e-10,...
    ConstraintTolerance=1e-10,...
    InitBarrierParam=0.1,...
    MaxIterations=500);

Aineq=[];
bineq=[];
Aeq=[];
beq=[];
lb=[];
ub=[];
%nonlcon=



[Aeq,beq]=AssembleLuvSSTREAM(CtrlVar,MUA,BCs) ;


tLSQ=tic;
[x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,Aineq,bineq,Aeq,beq,@mycon,options);
%[x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,Aineq,bineq,Aeq,beq);
%[x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
tLSQ=toc(tLSQ);

ub=x(1:MUA.Nnodes);
vb=x(MUA.Nnodes+1:end);


UaPlots(CtrlVar,MUA,F,"-uv-",FigureTitle="Ua uv solve")
UaPlots(CtrlVar,MUA,F,[ub vb],FigureTitle="LSQ uv solve")
UaPlots(CtrlVar,MUA,F,[ub-F.ub vb-F.vb],FigureTitle="diff")

%%
function [c,ceq,gradc,gradceq]=mycon(x)

c=0;
ceq=0;
gradc=sparse(length(x),1);
gradceq=sparse(length(x),1);



end
