

% n=200 ; A = sprandsym(n,.1,0.1);  x=ones(n,1) ;


% 

% delete(gcp('nocreate')); parpool('Threads',8)

ParPool = gcp('nocreate') ;

if isempty(ParPool)
    parpool('Processes',8)
end
parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:genericError');
parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:SaveNotSupported');
warning('off','MATLAB:decomposition:genericError')
 
Solving="-uvh-" ;

% Set output files directory
[~,hostname]=system('hostname') ;
if contains(hostname,"DESKTOP-G5TCRTD")  % office Dell

    UserVar.ResultsFileDirectory="F:\Runs\Calving\PIG-TWG\ResultsFiles\";
    UserVar.InverseRestartFileDirectory="F:\Runs\Calving\PIG-TWG\InverseRestartFiles\";
    UserVar.InversionFileDirectory="F:\Runs\Calving\PIG-TWG\InversionFiles\";
    UserVar.MeshFileDirectory="F:\Runs\Calving\PIG-TWG\MeshFiles\";
    UserVar.ForwardRestartFileDirectory="F:\Runs\Calving\PIG-TWG\RestartFiles\";

elseif contains(hostname,"DESKTOP-BU2IHIR")   % home

    UserVar.ResultsFileDirectory="D:\Runs\Calving\PIG-TWG\ResultsFiles\";
    UserVar.InverseRestartFileDirectory="D:\Runs\Calving\PIG-TWG\InverseRestartFiles\";
    UserVar.InversionFileDirectory="D:\Runs\Calving\PIG-TWG\InversionFiles\";
    UserVar.MeshFileDirectory="D:\Runs\Calving\PIG-TWG\MeshFiles\";
    UserVar.ForwardRestartFileDirectory="D:\Runs\Calving\PIG-TWG\RestartFiles\";
    
elseif contains(hostname,"C23000099")   % home

    UserVar.ResultsFileDirectory="E:\Runs\Calving\PIG-TWG\ResultsFiles\";
    UserVar.InverseRestartFileDirectory="E:\Runs\Calving\PIG-TWG\InverseRestartFiles\";
    UserVar.InversionFileDirectory="E:\Runs\Calving\PIG-TWG\InversionFiles\";
    UserVar.MeshFileDirectory="E:\Runs\Calving\PIG-TWG\MeshFiles\";
    UserVar.ForwardRestartFileDirectory="E:\Runs\Calving\PIG-TWG\RestartFiles\";
else
    UserVar.ResultsFileDirectory=pwd+"\ResultsFiles\";
end





%load(UserVar.InverseRestartFileDirectory+"InverseRestartFile-Joughin-Ca1-Cs100000-Aa1-As100000-5km-Alim-Clim-.mat","CtrlVarInRestartFile","RunInfo","MUA","F","BCs","l")

load(UserVar.InverseRestartFileDirectory+"InverseRestartFile-Cornford-Ca1-Cs100000-Aa1-As100000-5km-Alim-Clim-.mat","CtrlVarInRestartFile","RunInfo","MUA","F","BCs","l")

load(UserVar.ForwardRestartFileDirectory+"Restart-FT-P-Duvh-TWIS-MR4-SM-TM001-Cornford-2k5km-Alim-Clim-Ca1-Cs100000-Aa1-As100000-InvMR5","CtrlVarInRestartFile","RunInfo","MUA","F","BCs","l")




CtrlVar=CtrlVarInRestartFile;
CtrlVar.InfoLevelNonLinIt=1;  CtrlVar.InfoLevel=1;

CtrlVar.uvGroupAssembly=false; CtrlVar.uvhGroupAssembly=false; CtrlVar.etaZero=10;

CtrlVar.Parallel.uvAssembly.spmd.nWorkers=[];

CtrlVar.Parallel.uvAssembly.spmd.isOn=true;
CtrlVar.Parallel.uvAssembly.parfeval.isOn=false;
CtrlVar.Parallel.isTest=true;



MUAworkers=[]; 

if contains(Solving,"-uv-")
    % F.ub=F.ub*0 ; F.vb=F.vb*0;
    [UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,MUAworkers) ;

    UaPlots(CtrlVar,MUA,F,"-uv-")

end



MUA.workers=[]; 
MUA.workers=BuildMuaWorkers(CtrlVar,MUA,MUA.workers) ;

% CtrlVar.uvhMatrixAssembly.ZeroFields=false; CtrlVar.uvhMatrixAssembly.Ronly=false;
% [UserVar,RunInfo,R2,K2]=uvhAssemblySPMD2(UserVar,RunInfo,CtrlVar,MUA,F,F);


if contains(Solving,"-uvh-")
    %% uvh
    CtrlVar.dt=0.001;
    [UserVar,RunInfo,F1,l1,BCs1,dt]=uvh(UserVar,RunInfo,CtrlVar,MUA,F,F,l,l,BCs) ;

end