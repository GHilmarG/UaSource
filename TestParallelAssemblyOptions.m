


NumWorkers=8 ;

ParPool = gcp('nocreate') ;

if isempty(ParPool)

    parpool('Processes',NumWorkers)

elseif (ParPool.NumWorkers~=NumWorkers)

    delete(gcp('nocreate'))
    parpool('Processes',NumWorkers)

end

parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:genericError');
parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:SaveNotSupported');
warning('off','MATLAB:decomposition:genericError')
warning('off','MATLAB:decomposition:LoadNotSupported') 


Solving="-uv-" ;


UserVar=FileDirectories;


%load(UserVar.InverseRestartFileDirectory+"InverseRestartFile-Joughin-Ca1-Cs100000-Aa1-As100000-5km-Alim-Clim-.mat","CtrlVarInRestartFile","RunInfo","MUA","F","BCs","l")

if contains(Solving,"-uv-")
    load(UserVar.InverseRestartFileDirectory+"InverseRestartFile-Cornford-Ca1-Cs100000-Aa1-As100000-5km-Alim-Clim-.mat","CtrlVarInRestartFile","RunInfo","MUA","F","BCs","l")
else
    load(UserVar.ForwardRestartFileDirectory+"Restart-FT-P-Duvh-TWIS-MR4-SM-TM001-Cornford-2k5km-Alim-Clim-Ca1-Cs100000-Aa1-As100000-InvMR5","CtrlVarInRestartFile","RunInfo","MUA","F","BCs","l")
end

load(UserVar.ForwardRestartFileDirectory+"Restart-FT-P-Duvh-TWIS-MR4-SM-TM001-Cornford-2k5km-Alim-Clim-Ca1-Cs100000-Aa1-As100000-InvMR5","CtrlVarInRestartFile","RunInfo","MUA","F","BCs","l")

tic
MUA.dM=decomposition(MUA.M,'chol','upper') ;
toc

RunInfo=UaRunInfo;  % reset runinfo


CtrlVar=CtrlVarInRestartFile;
CtrlVar.InfoLevelNonLinIt=1;  CtrlVar.InfoLevel=1;

CtrlVar.uvGroupAssembly=false; CtrlVar.uvhGroupAssembly=false; CtrlVar.etaZero=10;

CtrlVar.Parallel.uvAssembly.spmd.nWorkers=[];


CtrlVar.Parallel.Options="-none-" ;
CtrlVar.Parallel.Options="-auto-" ;

if CtrlVar.Parallel.Options=="-auto-"

    CtrlVar.Parallel.uvAssembly.spmd.isOn=true;
    CtrlVar.Parallel.uvAssembly.parfeval.isOn=false;
    CtrlVar.Parallel.uvhAssembly.spmd.isOn=true;
    CtrlVar.Parallel.Distribute=true ;

elseif CtrlVar.Parallel.Options=="-none-"

    CtrlVar.Parallel.uvAssembly.spmd.isOn=false;
    CtrlVar.Parallel.uvAssembly.parfeval.isOn=false;
    CtrlVar.Parallel.uvhAssembly.spmd.isOn=false;
    CtrlVar.Parallel.Distribute=false ;
    
end

CtrlVar.Parallel.isTest=false;


MUA=UpdateMUA(CtrlVar,MUA) ;


%% modified NR options
CtrlVar.ModifiedNRuvIntervalCriterion=4  ; CtrlVar.ModifiedNRuvReductionCriterion=0.5 ;    % Total time=237.335  Solver=49.1201 	 Assembly=92.7946 
% CtrlVar.ModifiedNRuvIntervalCriterion=2  ; CtrlVar.ModifiedNRuvReductionCriterion=0.5 ;  % Total time=193.215  Solver=53.4797 	 Assembly=79.4635
% CtrlVar.ModifiedNRuvIntervalCriterion=1  ; CtrlVar.ModifiedNRuvReductionCriterion=0.5 ;  % Total time=202.01 	 Solver=71.846 	     Assembly=78.921 
 CtrlVar.InfoLevelNonLinIt=5 ;
%%

if contains(Solving,"-uv-")
    % F.ub=F.ub*0 ; F.vb=F.vb*0;
    tTotal=tic;
    [UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l) ;
    tTotal=toc(tTotal);

    fprintf("Total time=%g \t Solver=%g \t Assembly=%g \n",tTotal,RunInfo.CPU.Solution.uv,RunInfo.CPU.Assembly.uv)

    UaPlots(CtrlVar,MUA,F,"-uv-")


    % Total time=5.5462 	 Solver=2.02292 	 Assembly=3.44374   C23000099        12  SPMD   Distributed
    % Total time=6.41722 	 Solver=2.10235 	 Assembly=4.25054   C23000099         8  SPMD   Distributed

    % Total time=17.3182 	 Solver=6.08226 	 Assembly=11.1686   C23000099           ~SPMD  ~Distributed

end




if contains(Solving,"-uvh-")
    %% uvh
    CtrlVar.dt=0.001;
    
    tTotal=tic;
    [UserVar,RunInfo,F1,l1,BCs1,dt]=uvh(UserVar,RunInfo,CtrlVar,MUA,F,F,l,l,BCs) ;
    tTotal=toc(tTotal);

    fprintf("Total time=%g \t Solver=%g \t Assembly=%g \n",tTotal,RunInfo.CPU.Solution.uvh,RunInfo.CPU.Assembly.uvh)

    % Total time=89.8162 	 Solver=47.6165 	 Assembly=22.6194       C23000099        SPMD(24)

    % Total time=94.3706 	 Solver=49.6399 	 Assembly=23.7048       C23000099        SPMD(12)
    % Total time=75.4851 	 Solver=30.5031 	 Assembly=24.2612       C23000099        SPMD(12)    distributed

    % Total time=141.966 	 Solver=47.9826 	 Assembly=52.8219       C23000099        SPMD(4)
    % Total time=134.659 	 Solver=28.9537 	 Assembly=68.6005       C23000099       ~SPMD(8)     distributed

    % Total time=158.332 	 Solver=51.4621 	 Assembly=68.7052       C23000099       ~SPMD(12)   ~distributed




    % Total time=116.786 	 Solver=62.4066 	 Assembly=29.1553       DESKTOP-BU2IHIR  SPMD(12)
    % Total time=174.072 	 Solver=62.0067 	 Assembly=69.6168       DESKTOP-BU2IHIR ~SPMD(12)
    
     UaPlots(CtrlVar,MUA,F,"-uv-")
    

end