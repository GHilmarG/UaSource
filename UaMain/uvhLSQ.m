



function [UserVar,RunInfo,F1,l1,BCs1,dt]=uvhLSQ(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)


dt=CtrlVar.dt;
x=[F1.ub;F1.vb;F1.h];

MLC=BCs2MLC(CtrlVar,MUA,BCs1);
if numel(l1.ubvb)~=numel(MLC.ubvbRhs) ; l1.ubvb=zeros(numel(MLC.ubvbRhs),1) ; end
if numel(l1.h)~=numel(MLC.hRhs) ; l1.h=zeros(numel(MLC.hRhs),1) ; end


[Aeq,beq,l]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs1,l1);


if isempty(l)
    l=beq*0 ;
end

%% find a scale for the problem

CtrlVar.uvhMatrixAssembly.ZeroFields=true ; CtrlVar.uvhMatrixAssembly.Ronly=true ;
[UserVar,RunInfo,R0]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);
scale=sqrt(1/(R0'*R0/2));

%%

fun = @(x) uvhRK(x,UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1,scale) ;


Algorithm="-MATLAB-levenberg-marquardt";
%Algorithm="-MATLAB-interior-point-";
Algorithm="-UaLSQ-";
%Algorithm="-compare-" ;

if contains(Algorithm,"-UaLSQ-")

    

        CtrlVar.InfoLevelBackTrack=0;  CtrlVar.InfoLevelNonLinIt=1 ;  CtrlVar.doplots=1 ;
        [x,l,R2,r2,Qslope0,dxNorm,dlNorm,residual,g,h,outputUa] = lsqNewton(CtrlVar,fun,x,l,Aeq,beq) ;

elseif contains(Algorithm,"-MATLAB-")

    %options.PlotFcn={@optimplotfirstorderoptUa,@optimplotresnormUa};

    % on an unconstrained problem (RadialIceCap) the levenberg-marquard  and trust-region reflective work well, but it is
    % surprisingly difficult to get the interior-point to work as well
    %
    % options.Algorithm="interior-point" ;
    % options.Algorithm="trust-region-reflective";

    options = optimoptions('lsqnonlin','Display','iter','MaxIterations',25,'SpecifyObjectiveGradient',true);
    lb=[] ; ub=[] ; A=[] ; b=[] ;   nonlcon=[] ;
    options.OptimalityTolerance = 1.000000e-15 ;
    options.StepTolerance = 1.000000e-15 ;
    options.FunctionTolerance=1e-10;

    if contains(Algorithm,"-MATLAB-levenberg-marquardt")

        options.Algorithm="levenberg-marquardt"; options.ScaleProblem='jacobian'; options.InitDamping=1e-12;



    elseif contains(Algorithm,"-MATLAB-interior-point")

        options.Algorithm="interior-point";  
        options.BarrierParamUpdate="monotone" ;
        options.InitBarrierParam=0.0001;
        options.SubproblemAlgorithm="factorization" ;
        options.OptimalityTolerance = 1.000000e-12 ;

    end


    [x,resnorm,residual,exitflag,outputMAT] = lsqnonlin(fun,x,lb,ub,A,b,Aeq,beq,nonlcon,options);

elseif contains(Algorithm,"-compare-")

    CtrlVar.InfoLevelBackTrack=0;  CtrlVar.InfoLevelNonLinIt=1 ;  CtrlVar.doplots=1 ;
    tUa=tic;
    [xUa,l,R2,r2,Qslope0,dxNorm,dlNorm,residual,g,h,outputUa] = lsqNewton(CtrlVar,fun,x,l,Aeq,beq) ;
    tUa=toc(tUa);



    options = optimoptions('lsqnonlin','Display','none','MaxIterations',25,'SpecifyObjectiveGradient',true);
    options.OptimalityTolerance = 1.000000e-15 ;
    options.StepTolerance = 1.000000e-15 ;
    lb=[] ; ub=[] ; A=[] ; b=[] ;   nonlcon=[] ;


    options.Algorithm="levenberg-marquardt"; options.ScaleProblem='jacobian'; options.InitDamping=1e-12;
    tMATLM=tic;
    [xMATLM,resnormLM,residualLM,exitflag,outputMATLM] = lsqnonlin(fun,x,lb,ub,A,b,Aeq,beq,nonlcon,options);
    tMATLM=toc(tMATLM);

    options = optimoptions('lsqnonlin','algorithm','interior-point','Display','none','MaxIterations',25,'SpecifyObjectiveGradient',true);
    options.BarrierParamUpdate="monotone" ; options.InitBarrierParam=1e10; options.SubproblemAlgorithm="cg";
    options.ScaleProblem='none';
    %options.OptimalityTolerance = 1.000000e-12 ;
    tMATIP=tic;
    [xMATIP,resnormIP,residual,exitflag,outputMATIP] = lsqnonlin(fun,x,lb,ub,A,b,Aeq,beq,nonlcon,options);
    tMATIP=toc(tMATIP);

    fprintf("    Ua: CPU=%g \t Resnorm=%g \t 1st-Opt=%g \t it=%i \n",tUa,R2,r2,outputUa.nIt)
    fprintf("MAT-LM: CPU=%g \t Resnorm=%g \t 1st-Opt=%g \t it=%i \n",tMATLM,resnormLM,outputMATLM.firstorderopt,outputMATLM.iterations)
    fprintf("MAT-IP: CPU=%g \t Resnorm=%g \t 1st-Opt=%g \t it=%i \n",tMATIP,resnormIP,outputMATIP.firstorderopt,outputMATIP.iterations)

    x=xUa;


end

n=MUA.Nnodes;
F1.ub=x(1:n) ;
F1.vb=x(n+1:2*n) ;
F1.h=x(2*n+1:3*n) ;




end





