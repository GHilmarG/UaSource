

% n=200 ; A = sprandsym(n,.1,0.1);  x=ones(n,1) ;


% delete(gcp('nocreate')); parpool('Processes',8)

% delete(gcp('nocreate')); parpool('Threads',8)



% Set output files directory
[~,hostname]=system('hostname') ;
if contains(hostname,"DESKTOP-G5TCRTD")  % office Dell

    UserVar.ResultsFileDirectory="F:\Runs\Calving\PIG-TWG\ResultsFiles\";
    UserVar.InverseRestartFileDirectory="F:\Runs\Calving\PIG-TWG\InverseRestartFiles\";
    UserVar.InversionFileDirectory="F:\Runs\Calving\PIG-TWG\InversionFiles\";
    UserVar.MeshFileDirectory="F:\Runs\Calving\PIG-TWG\MeshFiles\";
    UserVar.ForwardRestartFileDirectory="F:\Runs\Calving\PIG-TWG\RestartFiles";

elseif contains(hostname,"DESKTOP-BU2IHIR")   % home

    UserVar.ResultsFileDirectory="D:\Runs\Calving\PIG-TWG\ResultsFiles\";
    UserVar.InverseRestartFileDirectory="D:\Runs\Calving\PIG-TWG\InverseRestartFiles\";
    UserVar.InversionFileDirectory="D:\Runs\Calving\PIG-TWG\InversionFiles\";
    UserVar.MeshFileDirectory="D:\Runs\Calving\PIG-TWG\MeshFiles\";
    UserVar.ForwardRestartFileDirectory="D:\Runs\Calving\PIG-TWG\RestartFiles";
    
elseif contains(hostname,"C23000099")   % home

    UserVar.ResultsFileDirectory="E:\Runs\Calving\PIG-TWG\ResultsFiles\";
    UserVar.InverseRestartFileDirectory="E:\Runs\Calving\PIG-TWG\InverseRestartFiles\";
    UserVar.InversionFileDirectory="E:\Runs\Calving\PIG-TWG\InversionFiles\";
    UserVar.MeshFileDirectory="E:\Runs\Calving\PIG-TWG\MeshFiles\";
    UserVar.ForwardRestartFileDirectory="E:\Runs\Calving\PIG-TWG\RestartFiles";
else
    UserVar.ResultsFileDirectory=pwd+"\ResultsFiles\";
end





%load(UserVar.InverseRestartFileDirectory+"InverseRestartFile-Joughin-Ca1-Cs100000-Aa1-As100000-5km-Alim-Clim-.mat","CtrlVarInRestartFile","RunInfo","MUA","F","BCs","l")

load(UserVar.InverseRestartFileDirectory+"InverseRestartFile-Cornford-Ca1-Cs100000-Aa1-As100000-5km-Alim-Clim-.mat","CtrlVarInRestartFile","RunInfo","MUA","F","BCs","l")

CtrlVar=CtrlVarInRestartFile;
CtrlVar.InfoLevelNonLinIt=1;  CtrlVar.InfoLevel=1;
CtrlVar.uvGroupAssembly=false;
CtrlVar.uvhGroupAssembly=false;

CtrlVar.Parallel.uvAssembly.spmd.nWorkers=[]; 

CtrlVar.Parallel.uvAssembly.spmd.isOn=false; 
CtrlVar.Parallel.uvAssembly.parfeval.isOn=false;
CtrlVar.Parallel.isTest=false; 


CtrlVar.etaZero=10; 

% F.ub=F.ub*0 ; F.vb=F.vb*0; 
[UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l) ; 

UaPlots(CtrlVar,MUA,F,"-uv-")

