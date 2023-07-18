function UserVar=Ua2D(UserVar,CtrlVarOnInput,varargin)

%% Driver for the 2HD Úa model
% 


if nargin==0
    UserVar=[]; 
    CtrlVarOnInput=[];
elseif nargin==1
    CtrlVarOnInput=[];
end



SetUaPath() %% 


warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
warning('off','MATLAB:decomposition:SaveNotSupported')
warning('off','MATLAB:decomposition:genericError')
 parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:genericError');
 parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:SaveNotSupported');

 
%% initialize some variables
RunInfo=UaRunInfo; 
Fm1=UaFields;
BCsAdjoint=BoundaryConditions;
Meas=Measurements;
Priors=PriorProbabilityDistribution;
InvStartValues=InversionValues;
InvFinalValues=InversionValues; 

RunInfo.Forward.AdaptiveTimeSteppingTimeStepModifiedForOutputs=0;
Lubvb=[];
Ruv=[];

WallTime0=tic;
%% Clear any persistent variables
ClearPersistentUaVariables();


%% Define CtrlVar

% The user can define the fields of CtrlVar by:
% 
% 1) Specifying the fields of CtlrVar in Ua2D_InitialUserInput.m
% 2) By handing over CtrlVar as a second argument to Ua
%
% If the user uses both of these options (the most typical situation) 
% then any fields defined on input to Ua take precedence.
%
% So for example, if Ua is called using the call
%
%   CtrlVar.dt=1; 
%   Ua(UserVar,CtrlVar)
%
% and within the Ua2D_InitialUserInput.m there is a line:
% 
%  CtrlVar.dt=10; 
%
% Then the value used for CtrlVar.dt is the one given as input to Ua, i.e. CtrlVar.dt=1
%
%


% get the Ua default values for the CtlrVar
CtrlVar=Ua2D_DefaultParameters();


if  ~isempty(CtrlVarOnInput)
    % now replace default CtrlVar fields with those of CtrlVarOnInput
    % of the input CtrlVarOnInput.
    fprintf("\n ===== The fields of CtrlVar given as input to Ua replace corresponding fields of CtrlVar defined in Ua2D_DefaultParameters.m \n")
    
    CtrlVar=ReplaceStructureFields(CtrlVar,CtrlVarOnInput);

end

[UserVar,CtrlVar]=GetInitialInputs(UserVar,CtrlVar,varargin{:});



if  ~isempty(CtrlVarOnInput)
    % now overwrite the changes to CtrlVar in Ua2D_InitialUserInput using the fields
    % of the input CtrlVarOnInput.
    fprintf("\n ===== The fields of CtrlVar given as input to Ua replace corresponding fields of CtrlVar defined in Ua2D_InitialUserInput.m \n")
    
    CtrlVar=ReplaceStructureFields(CtrlVar,CtrlVarOnInput);

end


%%

if ~isfield(CtrlVar,'fidlog')
    CtrlVar.fidlog=1;
end



% do some basic test on the vality of the CtrlVar fields, validate CtlrVar
CtrlVar=CtrlVarValidityCheck(CtrlVar);


[pathstr,name]=fileparts(CtrlVar.GmshFile); CtrlVar.GmshFile=[pathstr,name]; % get rid of eventual file extension



%%  Now initial information about the run has been defined
% write out some basic information about the type of run selected
PrintRunInfo(CtrlVar);
RunInfo.Message="Start of Run";
CtrlVar.RunInfoMessage=RunInfo.Message;
RunInfo.File.Name=CtrlVar.Experiment+"-RunInfo.txt";
if CtrlVar.WriteRunInfoFile
    if CtrlVar.Restart
        RunInfo.File.fid = fopen(RunInfo.File.Name,'a');
        fprintf(RunInfo.File.fid,'  Restart run starts on %s \n',datetime('now'));
    else
        RunInfo.File.fid = fopen(RunInfo.File.Name,'w');
    end
end

%% Get input data
if CtrlVar.InverseRun %  inverse run
    
    if CtrlVar.Restart %  inverse restart run
        
        
        RunInfo.Message="Getting inputs for an inverse restart run";
        CtrlVar.RunInfoMessage=RunInfo.Message;
        [UserVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=...
            GetInputsForInverseRestartRun(UserVar,CtrlVar,RunInfo);
        
    else % New inverse run
        
        
        RunInfo.Message="Getting inputs for a new inverse run";
        CtrlVar.RunInfoMessage=RunInfo.Message;
        CtrlVar.CurrentRunStepNumber=1; 
        % First get the usual input for a forward run
        
        [UserVar,RunInfo,MUA,BCs,F,l]=GetInputsForForwardRun(UserVar,CtrlVar,RunInfo);
        
        if CtrlVar.OnlyMeshDomainAndThenStop
            return
        end
        
        % now get the additional variables specific to an inverse run
        [UserVar,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=GetInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,RunInfo);
        
    end
    
else
    
    if CtrlVar.Restart %  forward restart run
        
        
        RunInfo.Message="Getting inputs for a forward restart run";
        CtrlVar.RunInfoMessage=RunInfo.Message;
        [UserVar,CtrlVarInRestartFile,MUA,BCs,F,l,RunInfo]=GetInputsForForwardRestartRun(UserVar,CtrlVar,RunInfo);
        
        
        % When reading the restart file the restart values of CtrlVar are all discarded,
        % however:
        CtrlVar.time=CtrlVarInRestartFile.time;   
        CtrlVar.RestartTime=CtrlVarInRestartFile.time;
       
        
        % Generally, I want dt used in the run to be the dt value in the restart file. However, it is possible that the restart file was
        % generated in a previous time independent run where dt was not defined, but the restart run is time dependent and the user has defined
        % dt in the DefineInitialInputs.m for that time dependent run. In that case I do want to use dt as defined in DefineInitialInputs.m
        if ~(~CtrlVarInRestartFile.TimeDependentRun && CtrlVar.TimeDependentRun)  
            % Generally the dt used should be the dt in the restart file. 
            CtrlVar.dt=CtrlVarInRestartFile.dt;
        end
        
        F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ; 
        CtrlVar.CurrentRunStepNumber=CtrlVarInRestartFile.CurrentRunStepNumber;
        
        clearvars time dt CurrentRunStepNumber
        
    else % New forward run (ie not a restart)
        
        
        RunInfo.Message="Getting inputs for a new forward run";
        CtrlVar.RunInfoMessage=RunInfo.Message;
        [UserVar,RunInfo,MUA,BCs,F,l]=GetInputsForForwardRun(UserVar,CtrlVar,RunInfo);
        
        if CtrlVar.OnlyMeshDomainAndThenStop
            return
        end
    end
    
end

MUA=UpdateMUA(CtrlVar,MUA); % Just in case something about the def of MUA has changed since creation of restart files
F.x=MUA.coordinates(:,1) ;  F.y=MUA.coordinates(:,2) ; F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ;  
%% RunInfo initialisation
RunInfo.Message="Start of Run";
CtrlVar.RunInfoMessage=RunInfo.Message;
RunInfo.File.Name=CtrlVar.Experiment+"-RunInfo.txt";

if CtrlVar.WriteRunInfoFile
    if CtrlVar.Restart
        if ~isnan(RunInfo.File.fid) ; fclose(RunInfo.File.fid) ; end
        RunInfo.File.fid = fopen(RunInfo.File.Name,'a');
        fprintf(RunInfo.File.fid,'  Restart run starts on %s \n',datetime('now'));
    else
        RunInfo.File.fid = fopen(RunInfo.File.Name,'w');
    end
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
RunInfo.Message="All initial inputs now defined.";  % this is a string, will only work correclty post Matlab 2017b.
CtrlVar.RunInfoMessage=RunInfo.Message;


F.h=F.s-F.b;
[F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow); 

% pointers to the elements of Boundary.Edges where u and v are fixed
% Boundary.uFixedEdgesAllPtrs=logical(prod(double(ismember(Boundary.Edges,ufixednode)')));
% Boundary.vFixedEdgesAllPtrs=logical(prod(double(ismember(Boundary.Edges,vfixednode)')));
% xint, yint :  matrices of coordinates of integration points. Nele  x nip
% Xint, Yint :  vectors of unique coordinates of integration points
% [DTxy,TRIxy,DTint,TRIint,Xint,Yint,xint,yint,Iint]=TriangulationNodesIntegrationPoints(MUA);


%%
if CtrlVar.doInverseStep   % -inverse
    
    
    RunInfo.Message="Start of inverse run.";
    CtrlVar.RunInfoMessage=RunInfo.Message;
     CtrlVar.CurrentRunStepNumber=1; 
    CtrlVar.DefineOutputsInfostring="Start of inverse run";
    CtrlVar.DefineOutputsCounter=1;
    InvFinalValues=InversionValues;
    fprintf(' Calling DefineOutputs. DefineOutputsInfostring=%s , DefineOutputsCounter=%i \n ',CtrlVar.DefineOutputsInfostring,CtrlVar.DefineOutputsCounter)
    UserVar=CreateOutputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
    
    
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
        InvertForModelParameters(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
    
    
    F.C=InvFinalValues.C          ; %fprintf(CtrlVar.fidlog,' C set equal to InvFinalValues.C. ');
    F.AGlen=InvFinalValues.AGlen  ; %fprintf(CtrlVar.fidlog,' AGlen set equal InvFinalValues.AGlen \n ');
    F.m=InvFinalValues.m ; 
    F.n=InvFinalValues.n ;
    
    [UserVar,RunInfo,F,l,drdu,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
    
    
    if CtrlVar.Inverse.WriteRestartFile
        
        WriteAdjointRestartFile(UserVar,CtrlVar,MUA,BCs,F,F.GF,l,RunInfo,InvStartValues,Priors,Meas,BCsAdjoint,InvFinalValues);
        
    end
    
    % PlotInverse    inverse plots
    if CtrlVar.doplots
        PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,l,F.GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
    end
    
    CtrlVar.DefineOutputsInfostring="End of Inverse Run";
    CtrlVar.DefineOutputsCounter=CtrlVar.DefineOutputsCounter+1;
    fprintf(' Calling DefineOutputs. DefineOutputsInfostring=%s , DefineOutputsCounter=%i \n ',CtrlVar.DefineOutputsInfostring,CtrlVar.DefineOutputsCounter)
    UserVar=CreateOutputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);

    SayGoodbye(CtrlVar,RunInfo);
    return  % This is the end of the (inverse) run
    
end

%% DefineOutputs
CtrlVar.DefineOutputsCounter=0;
if CtrlVar.CreateOutputsBeginningOfRun  && ~CtrlVar.Restart
    CtrlVar.DefineOutputsInfostring="First call";
    CtrlVar.DefineOutputsCounter=CtrlVar.DefineOutputsCounter+1;
    
    fprintf(' Calling DefineOutputs. DefineOutputsInfostring=%s , DefineOutputsCounter=%i \n ',CtrlVar.DefineOutputsInfostring,CtrlVar.DefineOutputsCounter)
    UserVar=CreateOutputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
    
    if CtrlVar.DefineOutputsCounter>=CtrlVar.DefineOutputsMaxNrOfCalls
        fprintf(' Exiting because number of calls to DefineOutputs (%i) >= CtrlVar.DefineOutputsMaxNrOfCalls (%i) /n',...
            CtrlVar.DefineOutputsCounter,CtrlVar.DefineOutputsMaxNrOfCalls)
        return
    end
end
%


CtrlVar.CurrentRunStepNumber0=CtrlVar.CurrentRunStepNumber;
CtrlVar.time0=CtrlVar.time;



RunInfo.Forward.IterationsTotal=0; 
RunInfo.Forward.uvhConverged=true; 
%%  RunStep Loop
while 1
    
    RunInfo.Message="-RunStepLoop-"; % While within run-step loop the Message field always contains the string "-RunStepLoop-"
    CtrlVar.RunInfoMessage=RunInfo.Message;
    RunInfo.CPU.WallTime=duration(0,0,toc(WallTime0));
    
    %% check run-step stop criteria
    if CtrlVar.CurrentRunStepNumber >=(CtrlVar.TotalNumberOfForwardRunSteps+CtrlVar.CurrentRunStepNumber0)
        
        fprintf('Exiting run-step loop because total number of steps reached. \n')
        %       fprintf('CtrlVar.CurrentRunStepNumber=%i \n',CtrlVar.CurrentRunStepNumber)
        %       fprintf('CtrlVar.CurrentRunStepNumber0=%i\n',CtrlVar.CurrentRunStepNumber0)
        break
    end
    
    if (CtrlVar.TotalTime - CtrlVar.time) <= CtrlVar.dtmin
        fprintf('Exiting time loop because total time reached. \n')
        break
    end
    
    if CtrlVar.TimeDependentRun && CtrlVar.dt < CtrlVar.dtmin  && ~RunInfo.Forward.AdaptiveTimeSteppingTimeStepModifiedForOutputs
        fprintf('Exiting run-step loop because time step too small (%g<%g)\n',CtrlVar.dt,CtrlVar.dtmin)
        TempFile=CtrlVar.Experiment+"-UaDumpTimeStepTooSmall.mat";
        fprintf(CtrlVar.fidlog,' saving variables in %s \n ',TempFile) ;
        save(TempFile,'-v7.3')
        break
    end
    
    
    
    if CtrlVar.UseUserDefinedRunStopCriterion
        
        [UserVar,isStop]=DefineRunStopCriterion(UserVar,RunInfo,CtrlVar,MUA,BCs,F) ;
        
        if isStop
            
            fprintf('Exiting run-step loop based on the user defined run stop criteria as specified in DefineRunStopCriteria.m \n')
            
            break
        end
        
    end
    
    %%
    

    CtrlVar.CurrentRunStepNumber=CtrlVar.CurrentRunStepNumber+1;
    

    if CtrlVar.InfoLevel >= 1 
        fprintf('\n \t ----------------------------------------> Current run step: %i <-------------------------------\n',CtrlVar.CurrentRunStepNumber) ;  
    end
    
    if CtrlVar.PlotWaitBar 
        multiWaitbar('Run steps','Value',(CtrlVar.CurrentRunStepNumber-1-CtrlVar.CurrentRunStepNumber0)/CtrlVar.TotalNumberOfForwardRunSteps);
        multiWaitbar('Time','Value',CtrlVar.time/CtrlVar.TotalTime);
    end
    
    MUA=UpdateMUA(CtrlVar,MUA);
    F.x=MUA.coordinates(:,1) ;  F.y=MUA.coordinates(:,2) ; 
    %% -adapt time step   automated time stepping 
    if CtrlVar.TimeDependentRun
        [RunInfo,CtrlVar.dt,CtrlVar.dtRatio]=AdaptiveTimeStepping(UserVar,RunInfo,CtrlVar,MUA,F);
    end
    
    
    %% [------------------adapt mesh    adaptive meshing,  adapt mesh, adapt-mesh
    if CtrlVar.AdaptMesh || CtrlVar.FEmeshAdvanceRetreat || CtrlVar.ManuallyDeactivateElements || CtrlVar.LevelSetMethodAutomaticallyDeactivateElements
        
        [UserVar,RunInfo,MUA,BCs,F,l]=AdaptMesh(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,Ruv,Lubvb);
        CtrlVar.AdaptMeshInitial=0;
        F.x=MUA.coordinates(:,1) ;  F.y=MUA.coordinates(:,2) ; 
        
        if MUA.Nele==0
            fprintf('FE mesh is empty \n ')
            break ;
        end
        
        if CtrlVar.doplots  && CtrlVar.PlotMesh
            figMesh=FindOrCreateFigure("Mesh");
            clf(figMesh) ; PlotMuaMesh(CtrlVar,MUA); hold on
            [xGL,yGL]=PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],'r','LineWidth',2);
            if ~isempty(xGL)
                Temp=figMesh.CurrentAxes.Title.String;
                figMesh.CurrentAxes.Title.String=[Temp(:)',{'Grounding line in red'}];
            end
            if ~isempty(F.LSF) && CtrlVar.LevelSetMethod    % Level Set  
                hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,'b','LineWidth',2) ;
                 Temp=figMesh.CurrentAxes.Title.String;
                figMesh.CurrentAxes.Title.String=[Temp(:)',{'Level set zero line in blue'}];
            end
            hold off
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
    
    %% [----------------------- Changes to F required at each RunStep external to the solver
    
    % There are questions as to which field should be updated manually first.
    % Geometry modifications might depend on the level-set, so call GetCalving before
    % GetGeometryAndDensities.
    
 
    [UserVar,F]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,F);
    [UserVar,F]=GetAGlenDistribution(UserVar,CtrlVar,MUA,F);
    
    if ~CtrlVar.doInverseStep
        if CtrlVar.TimeDependentRun
            [UserVar,F]=GetGeometryAndDensities(UserVar,CtrlVar,MUA,F,CtrlVar.GeometricalVarsDefinedEachTransienRunStepByDefineGeometry);
        else
            [UserVar,F]=GetGeometryAndDensities(UserVar,CtrlVar,MUA,F,CtrlVar.GeometricalVarsDefinedEachDiagnosticRunStepByDefineGeometry);
        end
    end

    if CtrlVar.UpdateBoundaryConditionsAtEachTimeStep
        [UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F);
        F=StartVelocity(CtrlVar,MUA,BCs,F);  % start velocity might be a function of GF
    end
 
    % get mass-balance after any modifications to geometry, as mass balance might depent
    % on geometry. 
    [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F);
    
    

    %%  -------------------------------------------------------------------------------------]

    %% "-uv-"
    %if ~CtrlVar.TimeDependentRun % Time independent run.  Solving for velocities for a given geometry (diagnostic step).
    if CtrlVar.UaRunType=="-uv-" % Time independent run.  Solving for velocities for a given geometry (diagnostic step).

        %% Diagnostic calculation (uv)
        if CtrlVar.InfoLevel >= 1 ; fprintf(CtrlVar.fidlog,' ==> Time independent step. Current run step: %i \n',CtrlVar.CurrentRunStepNumber) ;  end

        RunInfo.Message="-RunStepLoop- Diagnostic step. Solving for velocities.";
        CtrlVar.RunInfoMessage=RunInfo.Message;

        [UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);


    elseif CtrlVar.UaRunType=="-h-" % Time independent run.  Solving for velocities for a given geometry (diagnostic step).

        %% Diagnostic calculation (uv)
        if CtrlVar.InfoLevel >= 1 ; fprintf(CtrlVar.fidlog,' ==> Time independent step. Current run step: %i \n',CtrlVar.CurrentRunStepNumber) ;  end

        RunInfo.Message="-RunStepLoop- Diagnostic step. Solving for thickness.";
        CtrlVar.RunInfoMessage=RunInfo.Message;

        
        [UserVar,h,lambda]=hEquation(UserVar,CtrlVar,MUA,F,BCs);



        % "-uvh-" ; "-uv-h-"
    elseif CtrlVar.UaRunType=="-uvh-"  ||  CtrlVar.UaRunType=="-uv-h-"


        %        0  : values at t      This is F0
        %        1  : values at t+dt   This is F.
        %       at start, F is explicit guess for values at t+dt
        %       an end,   F are converged values at t+dt

        % RunInfo



        if numel(RunInfo.Forward.time) < CtrlVar.CurrentRunStepNumber 
            RunInfo.Forward.time=[RunInfo.Forward.time;RunInfo.Forward.time+NaN];
            RunInfo.Forward.dt=[RunInfo.Forward.dt;RunInfo.Forward.dt+NaN];
        end
        RunInfo.Forward.time(CtrlVar.CurrentRunStepNumber)=CtrlVar.time;
        RunInfo.Forward.dt(CtrlVar.CurrentRunStepNumber)=CtrlVar.dt;
        
        % "-uvh-"
        %
        % if CtrlVar.UaRunType=="-uvh-"  ||  CtrlVar.UaRunType=="-uv-h-"
        if CtrlVar.UaRunType=="-uvh-"
            
            
            RunInfo.Message="-RunStepLoop- Time dependent step. Solving implicitly for velocities and thickness.";
            CtrlVar.RunInfoMessage=RunInfo.Message;
            
            
            fprintf(...
                '\n ==========================>  Implicit uvh going from t=%-.10g to t=%-.10g with dt=%-g. Done %-g %% of total time, and  %-g %% of steps. (%s) \n ',...
                CtrlVar.time,CtrlVar.time+CtrlVar.dt,CtrlVar.dt,100*CtrlVar.time/CtrlVar.TotalTime,...
                100*(CtrlVar.CurrentRunStepNumber-1-CtrlVar.CurrentRunStepNumber0)/CtrlVar.TotalNumberOfForwardRunSteps,datetime('now'));
            
            if CtrlVar.WriteRunInfoFile
                
                RunInfo.CPU.Total=duration(0,0,cputime);
                RunInfo.CPU.WallTime=duration(0,0,toc(WallTime0));
                fprintf(RunInfo.File.fid,...
                    '  t-dt-tCPU-it-itTotal-date-uvhAssembly-uvhSolution-WallTime , %g , %g , %s  ,  %i , %i , %s , %s , %s , %s \n',...
                    CtrlVar.time,CtrlVar.dt,RunInfo.CPU.Total,RunInfo.Forward.uvhIterations(CtrlVar.CurrentRunStepNumber),...
                    RunInfo.Forward.IterationsTotal,datetime('now'),...
                    duration(0,0,RunInfo.CPU.Assembly.uvh),duration(0,0,RunInfo.CPU.Solution.uvh),RunInfo.CPU.WallTime);
            end
            
            
            if CtrlVar.InitialDiagnosticStep   % if not a restart step, and if not explicitly requested by user, then do not do an inital dignostic step
                %% diagnostic step, solving for uv.  Always needed at a start of a transient run. Also done if asked by the user.
                CtrlVar.InitialDiagnosticStep=0;
                
                fprintf(CtrlVar.fidlog,' initial diagnostic step at t=%-.15g \n ',CtrlVar.time);
                
                [UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
                
                
                %ub0=ub ; ud0=ud ; vb0=vb ; vd0=vd;
                
                
                if (ReminderFraction(CtrlVar.time,CtrlVar.DefineOutputsDt)<1e-5 || CtrlVar.DefineOutputsDt==0 )
                    CtrlVar.DefineOutputsInfostring="Diagnostic step";
                    CtrlVar.DefineOutputsCounter=CtrlVar.DefineOutputsCounter+1;
                    fprintf(' Calling DefineOutputs. DefineOutputsInfostring=%s , DefineOutputsCounter=%i \n ',CtrlVar.DefineOutputsInfostring,CtrlVar.DefineOutputsCounter)
                    
                    UserVar=CreateOutputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
                    if CtrlVar.DefineOutputsCounter>=CtrlVar.DefineOutputsMaxNrOfCalls
                        fprintf(' Exiting because number of calls to DefineOutputs (%i) >= CtrlVar.DefineOutputsMaxNrOfCalls (%i) /n',...
                            CtrlVar.DefineOutputsCounter,CtrlVar.DefineOutputsMaxNrOfCalls)
                        return
                    end
                end
            end
            
            % Now that the velocity has been calculated, we can ask for the calving parameters
           [UserVar,F]=GetCalving(UserVar,CtrlVar,MUA,F,BCs);  % Level Set  
            
            F0=F;  %
            
            
            %% get an explicit estimate for u, v and h at the end of the time step
            
            %
            % F0 is the converged solution from the previous time step
            % F0.dubdt is based on F0 and the previous solution to F0, which is referred to as Fm1, but is not saved
            % F0.dubdt=(F0.ub-Fm1.ub)/dt  (where dt is the time step between Fm1 and F0.)
            %
            
            [UserVar,RunInfo,F.ub,F.vb,F.ud,F.vd,F.h]=ExplicitEstimationForUaFields(UserVar,RunInfo,CtrlVar,MUA,F0,Fm1,BCs,l,BCs,l);

            
            %% advance the solution by dt using a fully implicit method with respect to u,v and h
            
            
            
            CtrlVar.time=CtrlVar.time+CtrlVar.dt;        % I here need the mass balance at the end of the time step, hence must increase t
            F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ; 
            [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F);
            CtrlVar.time=CtrlVar.time-CtrlVar.dt; % and then take it back to t at the beginning.
            F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ; 
            
            % uvh implicit step  (The F on input is based on an explicit estimate, on
            % return I have the implicit estimate. The explicit estimate is only there to
            % speed up the non-linear solver.
            [UserVar,RunInfo,F,l,BCs,dt]=uvh2(UserVar,RunInfo,CtrlVar,MUA,F0,F,l,l,BCs);
            
            CtrlVar.dt=dt;  % I might have changed dt within uvh
            
            if ~RunInfo.Forward.uvhConverged
                
                warning("Ua2D:WTSHTF","uvh did not converge")
                filename="DumpWTSHTD.mat";
                fprintf("Ua2D:Saving all variables in %s \n",filename)
                save(filename) 
                
                fprintf("Ua2D:calling WTSHTF\n")
                [UserVar,RunInfo,F,F0,l,Kuv,Ruv,Lubvb]= WTSHTF(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,Fm1,l);
                
            end
            

            CtrlVar.time=CtrlVar.time+CtrlVar.dt;
            F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ; 
            % Recalulating geometry based on floation not really needed here because uvh
            % does this implicitly.
            [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);
            [F,Fm1]=UpdateFtimeDerivatives(UserVar,RunInfo,CtrlVar,MUA,F,F0);
            
        % "-uv-h-"
        %elseif ~CtrlVar.Implicituvh % Semi-implicit time-dependent step. Implicit with respect to h, explicit with respect to u and v.
        elseif CtrlVar.UaRunType=="-uv-h-"
            
            RunInfo.Message="-RunStepLoop- Time dependent step. Solving explicitly for velocities and implicitly for thickness.";
            CtrlVar.RunInfoMessage=RunInfo.Message;
            

            fprintf(...
                '\n =========> Semi-Implicit uvh going from t=%-.10g to t=%-.10g with dt=%-g. Done %-g %% of total time, and  %-g %% of steps. (%s) \n ',...
                CtrlVar.time,CtrlVar.time+CtrlVar.dt,CtrlVar.dt,100*CtrlVar.time/CtrlVar.TotalTime,...
                100*(CtrlVar.CurrentRunStepNumber-1-CtrlVar.CurrentRunStepNumber0)/CtrlVar.TotalNumberOfForwardRunSteps,datetime('now'));

            tSemiImplicit=tic;                  % -uv


            F0=F;

            CtrlVar.time=CtrlVar.time+CtrlVar.dt;        % I here need the mass balance at the end of the time step, hence must increase t
            F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ;
            [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F);
            CtrlVar.time=CtrlVar.time-CtrlVar.dt; % and then take it back to t at the beginning.
            F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ;


            [UserVar,RunInfo,F,F0,l,Kuv,Ruv,Lubvb]= uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,Fm1,l);
            CtrlVar.InitialDiagnosticStep=0; 
            CtrlVar.time=CtrlVar.time+CtrlVar.dt; F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ;
            [F,Fm1]=UpdateFtimeDerivatives(UserVar,RunInfo,CtrlVar,MUA,F,F0);

            
            tSemiImplicit=toc(tSemiImplicit);
            if CtrlVar.InfoLevelCPU>=1
                fprintf(CtrlVar.fidlog,'SSTREAM semi-implicit step in %-g sec \n ',tSemiImplicit) ; 
            end
            
        end
        
        % update Level Set to current time using the new velocities
        if CtrlVar.LevelSetMethod
            [UserVar,RunInfo,F.LSF,F.LSFMask,F.LSFnodes,LSFlambda,F.LSFqx,F.LSFqy]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F);  % Level Set
        end
    end   % CtrlVar.TimeDependentRun
    
    %% calculations for this rund step are now done, only some plotting/writing issues do deal with


    
  
    %% plotting results
    
    % DefineOutputs
    
    if (ReminderFraction(CtrlVar.time,CtrlVar.DefineOutputsDt) < (CtrlVar.dt/(10*CtrlVar.DefineOutputsDt))) || CtrlVar.DefineOutputsDt==0 
        %(ReminderFraction(CtrlVar.time,CtrlVar.DefineOutputsDt)<1e-5 || CtrlVar.DefineOutputsDt==0 )
        CtrlVar.DefineOutputsInfostring="inside transient loop and inside run-step loop";
        CtrlVar.DefineOutputsCounter=CtrlVar.DefineOutputsCounter+1;
        
        if CtrlVar.MassBalanceGeometryFeedback>0
            CtrlVar.time=CtrlVar.time+CtrlVar.dt;  % I here need the mass balance at the end of the time step, hence must increase t
            F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ; 
            [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F);
            CtrlVar.time=CtrlVar.time-CtrlVar.dt; % and then take it back to t at the beginning. 
            F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ; 
            %[UserVar,as,ab,dasdh,dabdh]=GetMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time+CtrlVar.dt,s,b,h,S,B,rho,rhow,GF);
        end
        
        fprintf(' Calling DefineOutputs. DefineOutputsInfostring=%s , DefineOutputsCounter=%i \n ',CtrlVar.DefineOutputsInfostring,CtrlVar.DefineOutputsCounter)
        
        UserVar=CreateOutputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
        
        
        if CtrlVar.DefineOutputsCounter>=CtrlVar.DefineOutputsMaxNrOfCalls
            fprintf(' Exiting because number of calls to DefineOutputs (%i) >= CtrlVar.DefineOutputsMaxNrOfCalls (%i) /n',...
                CtrlVar.DefineOutputsCounter,CtrlVar.DefineOutputsMaxNrOfCalls)
            return
        end
    end
    
    if CtrlVar.WriteRestartFile==1 && mod(CtrlVar.CurrentRunStepNumber,CtrlVar.WriteRestartFileInterval)==0
        WriteForwardRunRestartFile(UserVar,CtrlVar,MUA,BCs,F,F.GF,l,RunInfo); 
    end
    
    
    
    
    
end

RunInfo.CPU.Total=duration(0,0,cputime);
RunInfo.Message="Calculations done. Creating outputs. ";
CtrlVar.RunInfoMessage=RunInfo.Message;

if CtrlVar.PlotWaitBar
    multiWaitbar('Run steps','Value',(CtrlVar.CurrentRunStepNumber-CtrlVar.CurrentRunStepNumber0)/(CtrlVar.TotalNumberOfForwardRunSteps-CtrlVar.CurrentRunStepNumber0));
    multiWaitbar('Time','Value',(CtrlVar.time-CtrlVar.time0) /(CtrlVar.TotalTime-CtrlVar.time0));
end


%% Possible final call to DefineOutputs



if CtrlVar.CreateOutputsEndOfRun
    CtrlVar.DefineOutputsInfostring="Last call";
    CtrlVar.DefineOutputsCounter=CtrlVar.DefineOutputsCounter+1;
    if CtrlVar.MassBalanceGeometryFeedback>0
        CtrlVar.time=CtrlVar.time+CtrlVar.dt; F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ; 
        [UserVar,F]=GetMassBalance(UserVar,CtrlVar,MUA,F);
        CtrlVar.time=CtrlVar.time-CtrlVar.dt; F.time=CtrlVar.time ;  F.dt=CtrlVar.dt ; 
        %[UserVar,as,ab,dasdh,dabdh]=GetMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time+CtrlVar.dt,s,b,h,S,B,rho,rhow,GF);
    end
    
    fprintf(' Calling DefineOutputs. DefineOutputsInfostring=%s , DefineOutputsCounter=%i \n ',CtrlVar.DefineOutputsInfostring,CtrlVar.DefineOutputsCounter)
    UserVar=CreateOutputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
    
    if CtrlVar.DefineOutputsCounter>=CtrlVar.DefineOutputsMaxNrOfCalls
        fprintf(' Exiting because number of calls to DefineOutputs (%i) >= CtrlVar.DefineOutputsMaxNrOfCalls (%i) /n',...
            CtrlVar.DefineOutputsCounter,CtrlVar.DefineOutputsMaxNrOfCalls)
        return
    end
end



%% saving outputs

if CtrlVar.WriteRestartFile==1 &&  mod(CtrlVar.CurrentRunStepNumber,CtrlVar.WriteRestartFileInterval)~=0
    WriteForwardRunRestartFile(UserVar,CtrlVar,MUA,BCs,F,F.GF,l,RunInfo); 
end

if CtrlVar.PlotWaitBar ;     multiWaitbar('CloseAll'); end

RunInfo.CPU.WallTime=duration(0,0,toc(WallTime0));
fprintf(CtrlVar.fidlog,' Wall-clock time : %s (hh:mm:ss) \n',RunInfo.CPU.WallTime) ;


if CtrlVar.fidlog~= 1 ; fclose(CtrlVar.fidlog); end


UserVar=DefineFinalReturnedValueOfUserVar(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);

SayGoodbye(CtrlVar,RunInfo)


end
