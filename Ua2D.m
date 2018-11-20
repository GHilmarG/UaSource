function Ua2D(UserVar,varargin)

%% Driver for the 2HD Úa model
% 



if nargin==0
    UserVar=[];
end

SetUaPath() %% set path

if ~exist(fullfile(cd,'Ua2D_InitialUserInput.m'),'file')
    
    fprintf('The input-file Ua2D_InitialUserInput.m not found in the working directory (%s).\n',pwd)
    fprintf('This input-file is required for Ua to run.\n')
    return
    
end

warning('off','MATLAB:triangulation:PtsNotInTriWarnId')

%% initialize some variables
RunInfo=UaRunInfo; 
Fm1=UaFields;
BCsAdjoint=BoundaryConditions;
Meas=Measurements;
Priors=PriorProbabilityDistribution;
InvStartValues=InversionValues;
InvFinalValues=InversionValues; 


Lubvb=[];
Ruv=[];

WallTime0=tic;
%% Clear any persistent variables
clear AdaptiveTimeStepping
clear AdaptMesh
clear BCs2MLC
clear CostFunctionValueAndGradient
clear EleAverageInterpolate
clear JGH
clear NrOfIntegrationPoints
clear LocalMeshRefinement
clear MeshAdvanceRetreat
clear Mesh2dEleSizeFunction
clear multiWaitbar
clear NewConjugatedGrad
% also those potentially defined in user input files
clear DefineSlipperyDistribution
clear DefineAGlenDistribution
clear DefineGeometry
clear DefineInputsForInverseRun
clear DefineDensities
clear DefineDesiredEleSize
clear DefineBoundaryConditions
clear DefineMassBalance
clear UaOutputs
% Ua utilities
clear PlotMeshScalarVariable 
clear PlotFEmesh

%% Define default values
CtrlVar=Ua2D_DefaultParameters();


%% Get user-defined parameter values
%  CtrlVar,UsrVar,Info,UaOuts
[UserVar,CtrlVar,MeshBoundaryCoordinates]=Ua2D_InitialUserInput(UserVar,CtrlVar,varargin{:});


%%
CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
clearvars MeshBoundaryCoordinates;


%%

if ~isfield(CtrlVar,'fidlog')
    CtrlVar.fidlog=1;
end



% do some basic test on the vality of the CtrlVar fields
CtrlVar=CtrlVarValidityCheck(CtrlVar);


[pathstr,name]=fileparts(CtrlVar.GmshFile); CtrlVar.GmshFile=[pathstr,name]; % get rid of eventual file extension



%%  Now initial information about the run has been defined
% write out some basic information about the type of run selected
PrintRunInfo(CtrlVar);
RunInfo.Message(1)="Start of Run";
RunInfo.File.Name=CtrlVar.Experiment+"-RunInfo.txt";
if CtrlVar.Restart
    RunInfo.File.fid = fopen(RunInfo.File.Name,'a');
    fprintf(RunInfo.File.fid,'  Restart run starts on %s \n',datetime('now'));
else
    RunInfo.File.fid = fopen(RunInfo.File.Name,'w');
end

%% Get input data
if CtrlVar.InverseRun %  inverse run
    
    if CtrlVar.Restart %  inverse restart run
        
        
        RunInfo.Message(numel(RunInfo.Message)+1)="Getting inputs for an inverse restart run";
        CtrlVar.RunInfoMessage=RunInfo.Message(end);
        [UserVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=...
            GetInputsForInverseRestartRun(UserVar,CtrlVar,RunInfo);
        
    else % New inverse run
        
        
        RunInfo.Message(numel(RunInfo.Message)+1)="Getting inputs for a new inverse run";
        CtrlVar.RunInfoMessage=RunInfo.Message(end);
        % First get the usual input for a forward run
        
        [UserVar,RunInfo,MUA,BCs,F,l,GF]=GetInputsForForwardRun(UserVar,CtrlVar,RunInfo);
        
        % now get the additional variables specific to an inverse run
        [UserVar,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=GetInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,GF,RunInfo);
        
    end
    
else
    
    if CtrlVar.Restart %  forward restart run
        
        
        RunInfo.Message(numel(RunInfo.Message)+1)="Getting inputs for a forward restart run";
        CtrlVar.RunInfoMessage=RunInfo.Message(end);
        [UserVar,CtrlVarInRestartFile,MUA,BCs,F,l,RunInfo]=GetInputsForForwardRestartRun(UserVar,CtrlVar,RunInfo);
        
        
        % When reading the restart file the restart values of CtrlVar are all discarded,
        % however:
        CtrlVar.time=CtrlVarInRestartFile.time;
        CtrlVar.RestartTime=CtrlVarInRestartFile.time;
        CtrlVar.dt=CtrlVarInRestartFile.dt;
        CtrlVar.CurrentRunStepNumber=CtrlVarInRestartFile.CurrentRunStepNumber;
        
        clearvars time dt CurrentRunStepNumber
        
    else % New forward run (ie not a restart)
        
        
        RunInfo.Message(numel(RunInfo.Message)+1)="Getting inputs for a new forward run";
        CtrlVar.RunInfoMessage=RunInfo.Message(end);
        [UserVar,RunInfo,MUA,BCs,F,l,GF]=GetInputsForForwardRun(UserVar,CtrlVar,RunInfo);
        
        if CtrlVar.OnlyMeshDomainAndThenStop
            return
        end
    end
    
end

if CtrlVar.TestUserInputs==1
    CtrlVar.TestUserInputs=0;
end


%% RunInfo initialisation
RunInfo.Message(1)="Start of Run";

RunInfo.File.Name=CtrlVar.Experiment+"-RunInfo.txt";

if CtrlVar.Restart
    if ~isnan(RunInfo.File.fid) ; fclose(RunInfo.File.fid) ; end
    RunInfo.File.fid = fopen(RunInfo.File.Name,'a');
    fprintf(RunInfo.File.fid,'  Restart run starts on %s \n',datetime('now'));
else
    RunInfo.File.fid = fopen(RunInfo.File.Name,'w');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    now all data specific to particular runs should have been defined %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  For convenience I assume that the user defines S, B, s and b.  The program
%  then calculates h=s-b and then s and b from h, B and S given the ice and
%  ocean specific density.  The thickness is preserved, and s and b are
%  consistent with the floating condition for a given ice tickness h, rho and
%  rhow.
if ~isfield(RunInfo,'Message') ; RunInfo.Message=[] ; end
RunInfo.Message(numel(RunInfo.Message)+1)="All initial inputs now defined.";  % this is a string, will only work correclty post Matlab 2017b.
CtrlVar.RunInfoMessage=RunInfo.Message(end);


F.h=F.s-F.b;
[F.b,F.s,F.h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);
%GF=GL2d(F.B,F.S,F.h,F.rhow,F.rho,MUA.connectivity,CtrlVar);



% pointers to the elements of Boundary.Edges where u and v are fixed
% Boundary.uFixedEdgesAllPtrs=logical(prod(double(ismember(Boundary.Edges,ufixednode)')));
% Boundary.vFixedEdgesAllPtrs=logical(prod(double(ismember(Boundary.Edges,vfixednode)')));
% xint, yint :  matrices of coordinates of integration points. Nele  x nip
% Xint, Yint :  vectors of unique coordinates of integration points
% [DTxy,TRIxy,DTint,TRIint,Xint,Yint,xint,yint,Iint]=TriangulationNodesIntegrationPoints(MUA);


%%
if CtrlVar.doInverseStep   % -inverse
    
    
    RunInfo.Message(numel(RunInfo.Message)+1)="Start of inverse run.";
    CtrlVar.RunInfoMessage=RunInfo.Message(end);
    
    CtrlVar.UaOutputsInfostring='Start of inverse run';
    CtrlVar.UaOutputsCounter=1;
    InvFinalValues=InversionValues;
    fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
    UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
    
    
    %         [db,dc] = Deblurr2D(NaN,...
    %             s,u,v,b,B,...
    %             sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,xMeas,yMeas,...
    %             Experiment,...
    %             coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,AGlen,C,...
    %             Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,nStep);
    %
    %x=coordinates(:,1); y=coordinates(:,2); DT = DelaunayTri(x,y); TRI=DT.Triangulation;
    %figure(21) ; trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,h) ;  title(' h')
        
    fprintf('\n =========================   Inverting for model parameters. =========================  \n')
    [UserVar,F,l,InvFinalValues,RunInfo]=...
        InvertForModelParameters(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
    
    
    F.C=InvFinalValues.C          ; %fprintf(CtrlVar.fidlog,' C set equal to InvFinalValues.C. ');
    F.AGlen=InvFinalValues.AGlen  ; %fprintf(CtrlVar.fidlog,' AGlen set equal InvFinalValues.AGlen \n ');
    F.m=InvFinalValues.m ; 
    F.n=InvFinalValues.n ;
    
    [UserVar,RunInfo,F,l,drdu,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
    
    if CtrlVar.doplots
         PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
    end
    
    if CtrlVar.Inverse.WriteRestartFile
        
        WriteAdjointRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo,InvStartValues,Priors,Meas,BCsAdjoint,InvFinalValues);
        
    end
    
    CtrlVar.UaOutputsInfostring='End of Inverse Run';
    CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
    fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
    UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);

    SayGoodbye(CtrlVar);
    return  % This is the end of the (inverse) run
    
end

%% UaOutputs
CtrlVar.UaOutputsCounter=0;
if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 )
    CtrlVar.UaOutputsInfostring='First call';
    CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
    
    fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
    UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
    
    if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls
        fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
            CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
        return
    end
end
%


CtrlVar.CurrentRunStepNumber0=CtrlVar.CurrentRunStepNumber;


RunInfo.Message(numel(RunInfo.Message)+1)="Run Step Loop.";
CtrlVar.RunInfoMessage=RunInfo.Message(end);
RunInfo.Forward.IterationsTotal=0; 
%%  RunStep Loop
while 1
    
    RunInfo.CPU.WallTime=duration(0,0,toc(WallTime0));
    
    if CtrlVar.CurrentRunStepNumber >=(CtrlVar.TotalNumberOfForwardRunSteps+CtrlVar.CurrentRunStepNumber0)
       
        fprintf('Exiting time loop because total number of steps reached. \n')
 %       fprintf('CtrlVar.CurrentRunStepNumber=%i \n',CtrlVar.CurrentRunStepNumber)
 %       fprintf('CtrlVar.CurrentRunStepNumber0=%i\n',CtrlVar.CurrentRunStepNumber0)
        break
    end
    
    if (CtrlVar.TotalTime - CtrlVar.time) <= CtrlVar.dtmin
        fprintf('Exiting time loop because total time reached. \n')
        break
    end
    
    if CtrlVar.TimeDependentRun && CtrlVar.dt <= CtrlVar.dtmin % I limit dt some small value for numerical reasons
        fprintf('Exiting time loop because time step too small (%g<%g)\n',CtrlVar.dt,CtrlVar.dtmin)
        TempFile=CtrlVar.Experiment+"-UaDumpTimeStepTooSmall.mat"; 
        fprintf(CtrlVar.fidlog,' saving variables in %s \n ',TempFile) ; 
        save(TempFile,'-v7.3')
        break
    end
    
    
    CtrlVar.CurrentRunStepNumber=CtrlVar.CurrentRunStepNumber+1;
    if CtrlVar.InfoLevel >= 1 
        fprintf('\n =========================================> Current run step: %i <==================================\n',CtrlVar.CurrentRunStepNumber) ;  
    end
    
    if CtrlVar.PlotWaitBar 
        multiWaitbar('Run steps','Value',(CtrlVar.CurrentRunStepNumber-1-CtrlVar.CurrentRunStepNumber0)/CtrlVar.TotalNumberOfForwardRunSteps);
        multiWaitbar('Time','Value',CtrlVar.time/CtrlVar.TotalTime);
    end
    
    MUA=UpdateMUA(CtrlVar,MUA);
    
    %% -adapt time step
    if CtrlVar.TimeDependentRun
        [RunInfo,CtrlVar.dt,CtrlVar.dtRatio]=AdaptiveTimeStepping(RunInfo,CtrlVar,CtrlVar.time,CtrlVar.dt);
    end
    
    
    %% [------------------adapt mesh    adaptive meshing,  adapt mesh, adapt-mesh
    if CtrlVar.AdaptMesh || CtrlVar.FEmeshAdvanceRetreat || CtrlVar.ManuallyDeactivateElements
        
        
        [UserVar,RunInfo,MUA,BCs,F,l,GF]=AdaptMesh(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,GF,Ruv,Lubvb);
        CtrlVar.AdaptMeshInitial=0;
        
        
        if MUA.Nele==0
            fprintf('FE mesh is empty \n ')
            break ;
        end
        
        if CtrlVar.doplots  && CtrlVar.PlotMesh
            
            FigMesh='Mesh';
            fig=findobj(0,'name',FigMesh);
            if isempty(fig)
                fig=figure('name',FigMesh);
                fig.Position=CtrlVar.PlotPosition;
            else
                fig=figure(fig);
                hold off
            end
            PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
        end

        
        
        if CtrlVar.AdaptMeshAndThenStop
            
            if ~isempty(CtrlVar.SaveAdaptMeshFileName)
                save(CtrlVar.SaveAdaptMeshFileName,'MUA') ;
                fprintf(' MUA was saved in %s .\n',CtrlVar.SaveAdaptMeshFileName);
            end
            fprintf('Exiting after remeshing because CtrlVar.AdaptMeshAndThenStop set to true. \n ')
            return
        end
        
    end
    
    %% -------------------------------------------------------------------------------------------]
    
    %% [----------------------- Changes to F required at each RunStep
    
    if ~CtrlVar.doInverseStep
        if CtrlVar.TimeDependentRun
            [UserVar,F,GF]=GetGeometryAndDensities(UserVar,CtrlVar,MUA,F,CtrlVar.GeometricalVarsDefinedEachTransienRunStepByDefineGeometry);
        else
            [UserVar,F,GF]=GetGeometryAndDensities(UserVar,CtrlVar,MUA,F,CtrlVar.GeometricalVarsDefinedEachDiagnosticRunStepByDefineGeometry);
        end
    end
    
    

    [UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F,GF);
    [UserVar,F]=GetAGlenDistribution(UserVar,CtrlVar,MUA,F,GF);

    if CtrlVar.UpdateBoundaryConditionsAtEachTimeStep
        [UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F,GF);
        F=StartVelocity(CtrlVar,MUA,BCs,F);  % start velocity might be a function of GF
    end
 
    
    [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F,GF);
    
    %%  -------------------------------------------------------------------------------------]
    
    
    if ~CtrlVar.TimeDependentRun % Time independent run.  Solving for velocities for a given geometry (diagnostic steo).
        
        %% Diagnostic calculation (uv)
        if CtrlVar.InfoLevel >= 1 ; fprintf(CtrlVar.fidlog,' ==> Time independent step. Current run step: %i \n',CtrlVar.CurrentRunStepNumber) ;  end
        
        RunInfo.Message(numel(RunInfo.Message)+1)="Diagnostic step. Solving for velocities.";
        CtrlVar.RunInfoMessage=RunInfo.Message(end);
        
        [UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
        
        
    else   % Time-dependent run
        
        %        0  : values at t      This is F0
        %        1  : values at t+dt   This is F.
        %       at start, F is explicit guess for values at t+dt
        %       an end,   F are converged values at t+dt
        
        if CtrlVar.Implicituvh % Fully implicit time-dependent step (uvh)
            
            
            RunInfo.Message(numel(RunInfo.Message)+1)="Time dependent step. Solving implicitly for velocities and thickness.";
            CtrlVar.RunInfoMessage=RunInfo.Message(end);
            
            
            fprintf(...
                '\n ---------> Implicit uvh going from t=%-.10g to t=%-.10g with dt=%-g. Done %-g %% of total time, and  %-g %% of steps. (%s) \n ',...
                CtrlVar.time,CtrlVar.time+CtrlVar.dt,CtrlVar.dt,100*CtrlVar.time/CtrlVar.TotalTime,100*(CtrlVar.CurrentRunStepNumber-1-CtrlVar.CurrentRunStepNumber0)/CtrlVar.TotalNumberOfForwardRunSteps,datetime('now'));
            
            if CtrlVar.WriteRunInfoFile
                
                RunInfo.CPU.Total=duration(0,0,cputime);
                RunInfo.CPU.WallTime=duration(0,0,toc(WallTime0));
                
                fprintf(RunInfo.File.fid,...
                    '  t-dt-tCPU-it-itTotal-date-uvhAssembly-uvhSolution-WallTime , %g , %g , %s  ,  %i , %i , %s , %s , %s , %s \n',...
                    CtrlVar.time,CtrlVar.dt,RunInfo.CPU.Total,RunInfo.Forward.Iterations,RunInfo.Forward.IterationsTotal,datetime('now'),...
                    duration(0,0,RunInfo.CPU.Assembly.uvh),duration(0,0,RunInfo.CPU.Solution.uvh),RunInfo.CPU.WallTime);
                
            end
            
            
            if CtrlVar.InitialDiagnosticStep   % if not a restart step, and if not explicitly requested by user, then do not do an inital dignostic step
                %% diagnostic step, solving for uv.  Always needed at a start of a transient run. Also done if asked by the user.
                CtrlVar.InitialDiagnosticStep=0;
                
                fprintf(CtrlVar.fidlog,' initial diagnostic step at t=%-.15g \n ',CtrlVar.time);
                
                [UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
                
                
                %ub0=ub ; ud0=ud ; vb0=vb ; vd0=vd;
                
                
                if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 )
                    CtrlVar.UaOutputsInfostring='Diagnostic step';
                    CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
                    fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
                    
                    UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
                    if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls
                        fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
                            CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
                        return
                    end
                end
            end
            
            F0=F;  % 
            
            
            %% get an explicit estimate for u, v and h at the end of the time step
            
            %
            % F0 is the converged solution from the previous time step
            % F0.dubdt is based on F0 and the previous solution to F0, which is referred to as Fm1, but is not saved 
            % F0.dubdt=(F0.ub-Fm1.ub)/dt  (where dt is the time step between Fm1 and F0.)
            %
            
            F=ExplicitEstimationForUaFields(CtrlVar,F,F0,Fm1);
            
            
            %% advance the solution by dt using a fully implicit method with respect to u,v and h
            uvhStep=1;
            while uvhStep==1  && CtrlVar.dt > CtrlVar.dtmin  % if uvh step does not converge, it is repeated with a smaller dt value
                CtrlVar.time=CtrlVar.time+CtrlVar.dt;        % I here need the mass balance at the end of the time step, hence must increase t
                [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F,GF);
                CtrlVar.time=CtrlVar.time-CtrlVar.dt; % and then take it back to t at the beginning.
                
                
                %Fguessed=F;  could use norm of difference between explicit and implicit to
                %control dt
                [UserVar,RunInfo,F,l,BCs,GF,dt]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F,l,l,BCs);
                
                CtrlVar.dt=dt;  % I might have changed dt within uvh
                
                if ~RunInfo.Forward.Converged
                    
                    uvhStep=1;  % continue within while loop
                    
                    filename="Dumpfile_Ua2D-"+CtrlVar.Experiment+".mat";
                    fprintf(' ===>>> uvh did not converge! Saving all data in a dumpfile %s \n',filename)
                    try
                        save(filename)
                    catch
                        warning('Ua2D:FileNotSaved',...
                            'could not save file %s.',filename)
                    end
                    
                    
                    fprintf(CtrlVar.fidlog,' =====>>> Reducing time step from %-g to %-g \n',CtrlVar.dt,CtrlVar.dt/10);
                    fprintf(CtrlVar.fidlog,'             Also resetting field and Lagrange variables. \n');
                    fprintf(CtrlVar.fidlog,'             Starting values for velocities at end of time step set by StartVelocity.m \n');
                    fprintf(CtrlVar.fidlog,'             Starting for s, b and h at end of time step set equal to values at beginning of time step. \n');
                    
                    CtrlVar.dt=CtrlVar.dt/10;
                    F.s=F0.s ; F.b=F0.b ; F.h=F0.h;
                    l.ubvb=l.ubvb*0; l.h=l.h*0;
                    F=StartVelocity(CtrlVar,MUA,BCs,F);
                    
                else
                    uvhStep=0;
                end
                
            end
            
            CtrlVar.time=CtrlVar.time+CtrlVar.dt;
            %CtrlVar.time=round(CtrlVar.time,14,'significant');
            
            [F.b,F.s,F.h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);  % This should not be needed as uvh already takes care of this.
           
            
            Fm1.dhdt=F0.dhdt ;
            Fm1.dubdt=F0.dubdt ; Fm1.dvbdt=F0.dvbdt;
            Fm1.duddt=F0.duddt ; Fm1.dvddt=F0.dvddt;
            
            if CtrlVar.dt==0
                F.dhdt=[];
                F.dubdt=[]; F.dvbdt=[];
                F.dsdt=[] ; F.dbdt=[];
            else
                F.dhdt=(F.h-F0.h)/CtrlVar.dt;
                F.dsdt=(F.s-F0.s)/CtrlVar.dt;
                F.dbdt=(F.b-F0.b)/CtrlVar.dt;
                F.dubdt=(F.ub-F0.ub)/CtrlVar.dt ; F.dvbdt=(F.vb-F0.vb)/CtrlVar.dt;
                F.duddt=(F.ud-F0.ud)/CtrlVar.dt ; F.dvddt=(F.vd-F0.vd)/CtrlVar.dt;
            end
            

            % at the beginning of next times step update:  F0=F1   
            
            
            
        elseif ~CtrlVar.Implicituvh % Semi-implicit time-dependent step. Implicit with respect to h, explicit with respect to u and v.
            
            
            RunInfo.Message(numel(RunInfo.Message)+1)="Time dependent step. Solving explicitly for velocities and implicitly for thickness.";
            CtrlVar.RunInfoMessage=RunInfo.Message(end);
            
            if CtrlVar.InfoLevel>0 ; fprintf(CtrlVar.fidlog,'Semi-implicit transient step. Advancing time from t=%-g to t=%-g \n',CtrlVar.time,CtrlVar.time+CtrlVar.dt);end
            
            %% Diagnostic calculation (uv)
            if CtrlVar.InfoLevel >= 1 ; fprintf(CtrlVar.fidlog,' ==> Diagnostic step (uv). Current run step: %i \n',CtrlVar.CurrentRunStepNumber) ;  end
            tSemiImplicit=tic;                  % -uv
            
            % Solving the momentum equations for uv at the beginning of the time interval
             
            F0=F;
            
            [UserVar,RunInfo,F,F0,l,Kuv,Ruv,Lubvb,GF]= uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,Fm1,l);

            Fm1.dhdt=F0.dhdt ;
            Fm1.dubdt=F0.dubdt ; Fm1.dvbdt=F0.dvbdt;
            Fm1.duddt=F0.duddt ; Fm1.dvddt=F0.dvddt;
                      
            CtrlVar.time=CtrlVar.time+CtrlVar.dt;
            CtrlVar.time=round(CtrlVar.time,14,'significant');
            
            %CtrlVar.hChange=1; CtrlVar.s=1; % h and s just changed
            tSemiImplicit=toc(tSemiImplicit);
            if CtrlVar.InfoLevel >= 1 && fprintf(CtrlVar.fidlog,'SSTREAM semi-implicit step in %-g sec \n ',tSemiImplicit) ; end
        end
    end
    
    %% calculations for this step are now done, only some plotting/writing issues do deal with
    
    % calculating derived quantities
    
    %GF = GL2d(F.B,F.S,F.h,F.rhow,F.rho,MUA.connectivity,CtrlVar);
    
  
    %% plotting results
    
    % UaOutputs
    
    if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 )
        CtrlVar.UaOutputsInfostring='inside transient loop and inside run-step loop';
        CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
        
        if CtrlVar.MassBalanceGeometryFeedback>0
            CtrlVar.time=CtrlVar.time+CtrlVar.dt;  % I here need the mass balance at the end of the time step, hence must increase t
            [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F,GF);
            CtrlVar.time=CtrlVar.time-CtrlVar.dt; % and then take it back to t at the beginning. 
            %[UserVar,as,ab,dasdh,dabdh]=GetMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time+CtrlVar.dt,s,b,h,S,B,rho,rhow,GF);
        end
        
        fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
        
        UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
        %UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,s,b,S,B,h,ub,vb,ud,vd,uo,vo,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,dasdh,dabdh,GF,BCs,l);
        
        if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls
            fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
                CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
            return
        end
    end
    
    if CtrlVar.WriteRestartFile==1 && mod(CtrlVar.CurrentRunStepNumber,CtrlVar.WriteRestartFileInterval)==0
        WriteForwardRunRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo); 
    end
    
    
    
    
    
end

RunInfo.CPU.Total=duration(0,0,cputime);
RunInfo.Message(numel(RunInfo.Message)+1)="Calculations done. Creating outputs. ";
CtrlVar.RunInfoMessage=RunInfo.Message(end);

if CtrlVar.PlotWaitBar
    multiWaitbar('Run steps','Value',(CtrlVar.CurrentRunStepNumber-CtrlVar.CurrentRunStepNumber0)/CtrlVar.TotalNumberOfForwardRunSteps);
    multiWaitbar('Time','Value',CtrlVar.time/CtrlVar.TotalTime);
end


%% plotting results

%[etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n);
%[wSurf,wSurfInt,wBedInt,wBed]=calcVerticalSurfaceVelocity(rho,rhow,h,S,B,b,ub,vb,as,ab,exx,eyy,xint,yint,MUA.coordinates,MUA.connectivity,MUA.nip,CtrlVar);

if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 )
    CtrlVar.UaOutputsInfostring='Last call';
    CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
    if CtrlVar.MassBalanceGeometryFeedback>0
        CtrlVar.time=CtrlVar.time+CtrlVar.dt;
        [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F,GF);
        CtrlVar.time=CtrlVar.time-CtrlVar.dt;
        %[UserVar,as,ab,dasdh,dabdh]=GetMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time+CtrlVar.dt,s,b,h,S,B,rho,rhow,GF);
    end
    
    fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
    UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
    %UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,s,b,S,B,h,ub,vb,ud,vd,uo,vo,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,dasdh,dabdh,GF,BCs,l);
    if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls
        fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
            CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
        return
    end
end



%% saving outputs

if CtrlVar.WriteRestartFile==1 &&  mod(CtrlVar.CurrentRunStepNumber,CtrlVar.WriteRestartFileInterval)~=0
    WriteForwardRunRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo); 
end

if CtrlVar.PlotWaitBar ;     multiWaitbar('CloseAll'); end

RunInfo.CPU.WallTime=duration(0,0,toc(WallTime0));
fprintf(CtrlVar.fidlog,' Wall-clock time : %s (hh:mm:ss) \n',RunInfo.CPU.WallTime) ;


if CtrlVar.fidlog~= 1 ; fclose(CtrlVar.fidlog); end

fclose(RunInfo.File.fid);

SayGoodbye(CtrlVar)


end
