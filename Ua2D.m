function Ua2D(UserRunParameters)

%% Driver for the 2HD Úa model
% Ua2D(UserRunParameters)
%
%

if nargin==0
    UserRunParameters=[];
end


SetUaPath() %% set path

if ~exist(fullfile(cd,'Ua2D_InitialUserInput.m'),'file')
    
    fprintf('The input-file Ua2D_InitialUserInput.m not found in the working directory (%s).\n',pwd)
    fprintf('This input-file is required for Ua to run.\n')
    return
    
end

warning('off','MATLAB:triangulation:PtsNotInTriWarnId')

%% initialize some variables
Info=UaRunInfo;
l=UaLagrangeVariables; % these are the Lagrange variables assosiated with boundary conditions.
% not to be confused with the Lagrange variables assosiated with solving the
% adjoint problem

RunInfo=[];
Lubvb=[];
GLdescriptors=[] ;
da0dt=[];
dsdt=NaN; dbdt=NaN; dhdt=NaN; Ruv=[];
dGFdt=[];  % get rid of this at a later stage
tTime=tic;
as=[];ab=[]; dasdh=[] ; dabdh=[];   BCs=[];
Experiment='';

%% Define default values
CtrlVar=Ua2D_DefaultParameters();
CtrlVar.UserParameters=UserRunParameters;


%% Get user-defined parameter values
%  CtrlVar,UsrVar,Info,UaOuts
[UserVar,CtrlVar,time,dt,MeshBoundaryCoordinates]=Ua2D_InitialUserInput(CtrlVar);

if ischar(UserVar)
    Experiment=UserVar ;
    CtrlVar.Experiment=Experiment;
end



%% copy Experiment, time, dt and MeshBoundaryCoordinates into CtrlVar
% and once that is done get rid of those.

CtrlVar.time=time;
CtrlVar.dt=dt;
CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
clearvars Experiment time dt MeshBoundaryCoordinates;

%%

if ~isfield(CtrlVar,'fidlog')
    CtrlVar.fidlog=1;
end

CtrlVar.MeshChanged=0;  % true if mesh changed in last adapt-meshing stage

% do some basic test on the vality of the CtrlVar fields
CtrlVar=CtrlVarValidityCheck(CtrlVar);

% if a log file is created, the this will be the name of the logfile
CtrlVar.Logfile=[CtrlVar.Experiment,'.log'];
[status,message]=copyfile(CtrlVar.Logfile,[CtrlVar.Logfile,'~']);  %  copy potential previous logfile

[pathstr,name]=fileparts(CtrlVar.GmshFile); CtrlVar.GmshFile=[pathstr,name]; % get rid of eventual file extension

%%  A RunInfo file is created and some basic information about the run is added to that file
%   Inspecting the contents of this run info file is sometimes a convenient way of
%   checking the status of the run.
if CtrlVar.Restart
    CtrlVar.InfoFile = fopen([CtrlVar.Experiment,'-RunInfo.txt'],'a');
else
    CtrlVar.InfoFile = fopen([CtrlVar.Experiment,'-RunInfo.txt'],'w');
end


%%  Now initial information about the run has been defined
% write out some basic information about the type of run selected
PrintRunInfo(CtrlVar);


%% Get input data
if ~CtrlVar.InverseRun %  forward run
    
    if CtrlVar.Restart  % Forward restart run
        
        [UserVar,MUA,BCs,time,dt,CurrentRunStepNumber,s,b,S,B,ub,vb,ud,vd,l,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,dasdh,dabdh,...
            dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1,GLdescriptors]=...
            GetInputsForForwardRestartRun(UserVar,CtrlVar);
        
        % When reading the restart file the restart values of CtrlVar are all discarded,
        % however:
        CtrlVar.time=time;
        CtrlVar.RestartTime=time;
        CtrlVar.dt=dt;
        CtrlVar.CurrentRunStepNumber=CurrentRunStepNumber;
        
        
        clearvars time dt CurrentRunStepNumber
        
        
    else % New forward run (ie not a restart)
        [UserVar,MeshChanged,MUA,BCs,s,b,S,B,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,dasdh,dabdh,...
            dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1]=...
            GetInputsForForwardRun(UserVar,CtrlVar);
        
        
        CtrlVar.MeshChanged=MeshChanged;
        
        if CtrlVar.OnlyMeshDomainAndThenStop
            return
        end
    end
    
else % inverse run
    
    if CtrlVar.Restart %  inverse restart run
        
        [UserVar,MUA,BCs,s,b,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,InvStartValues,Priors,Meas,BCsAdjoint,Info]=...
            GetInputsForInverseRestartRun(UserVar,CtrlVar);
        
    else % New inverse run
        
        [UserVar,MeshChanged,MUA,BCs,s,b,S,B,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,...
            dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1,...
            GF]=GetInputsForForwardRun(UserVar,CtrlVar);
        CtrlVar.MeshChanged=MeshChanged;
        
        
        
        % now get the additional variables specific to an inverse run
        [UserVar,InvStartValues,Priors,Meas,BCsAdjoint]=GetInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,CtrlVar.time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF);
        
    end
end

if CtrlVar.TestUserInputs==1
    CtrlVar.TestUserInputs=0;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    now all data specific to particular runs should have been defined %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  For convenience I assume that the user defines S, B, s and b.  The
%  program then calculates h=s-b and then s and b from h, B and S given the
%  ice and ocean specific density.  The thickness is preserved, and s and b
%  are consistent with the floating condition for a given ice tickness h, rho and rhow.

h=s-b;
[b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);



% pointers to the elements of Boundary.Edges where u and v are fixed
% Boundary.uFixedEdgesAllPtrs=logical(prod(double(ismember(Boundary.Edges,ufixednode)')));
% Boundary.vFixedEdgesAllPtrs=logical(prod(double(ismember(Boundary.Edges,vfixednode)')));
% xint, yint :  matrices of coordinates of integration points. Nele  x nip
% Xint, Yint :  vectors of unique coordinates of integration points
% [DTxy,TRIxy,DTint,TRIint,Xint,Yint,xint,yint,Iint]=TriangulationNodesIntegrationPoints(MUA);


%%
if CtrlVar.doInverseStep   % -inverse
    
    %         [db,dc] = Deblurr2D(NaN,...
    %             s,u,v,b,B,...
    %             sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,xMeas,yMeas,...
    %             Experiment,...
    %             coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,AGlen,C,...
    %             Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,nStep);
    %
    %x=coordinates(:,1); y=coordinates(:,2); DT = DelaunayTri(x,y); TRI=DT.Triangulation;
    %figure(21) ; trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,h) ;  title(' h')
    
    [UserVar,InvFinalValues,ub,vb,ud,vd,l,xAdjoint,yAdjoint,Info]=...
        InvertForModelParameters(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info);
    
    
    C=InvFinalValues.C          ; fprintf(CtrlVar.fidlog,' C set equal to InvFinalValues.C \n ');
    AGlen=InvFinalValues.AGlen  ; fprintf(CtrlVar.fidlog,' AGlen set equal InvFinalValues.AGlen \n ');
    m=InvFinalValues.m ; n=InvFinalValues.n ;
    
    
    % this calculation not really needed as AdjointNR2D should return converged ub,vb,ud,vd values for Cest and AGlenEst
    [UserVar,ub,vb,ud,vd,l,Kuv,Ruv,RunInfo,L]= uv(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF);
    
    
    if CtrlVar.doplots
        AdjointResultsPlots(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,...
            InvStartValues,Priors,Meas,BCsAdjoint,Info,InvFinalValues,xAdjoint,yAdjoint);
    end
    
    if CtrlVar.AdjointWriteRestartFile
        
        WriteAdjointRestartFile();
        
    end
    
    CtrlVar.UaOutputsInfostring='End of Inverse Run';
    CtrlVar.UaOutputsCounter=1;
    
    fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
    UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,s,b,S,B,h,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,dasdh,dabdh,GF,BCs,l);
    
    SayGoodbye(CtrlVar);
    return  % This is the end of the (inverse) run
    
end

%% UaOutputs
CtrlVar.UaOutputsCounter=0;
if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 )
    CtrlVar.UaOutputsInfostring='First call';
    CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
    
    fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
    
    UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,s,b,S,B,h,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,dasdh,dabdh,GF,BCs,l);
    if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls
        fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
            CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
        return
    end
end
%


CtrlVar.CurrentRunStepNumber0=CtrlVar.CurrentRunStepNumber;


%%  time loop
while 1
    
    if CtrlVar.CurrentRunStepNumber >=( CtrlVar.TotalNumberOfForwardRunSteps+CtrlVar.CurrentRunStepNumber0)
        fprintf('Exiting time loop because total number of steps reached. \n')
        break
    end
    
    if (CtrlVar.TotalTime - CtrlVar.time) <= CtrlVar.dtmin
        fprintf('Exiting time loop because total time reached. \n')
        break
    end
    
    if CtrlVar.TimeDependentRun && CtrlVar.dt <= CtrlVar.dtmin % I limit dt some small value for numerical reasons
        fprintf('Exiting time loop because time step too small (%g<%g)\n',CtrlVar.dt,CtrlVar.dtmin)
        TempFile=[CtrlVar.Experiment,'-UaDumpTimeStepTooSmall.mat']; fprintf(CtrlVar.fidlog,' saving variables in %s \n ',TempFile) ; save(TempFile,'-v7.3')
        break
    end
    
    
    CtrlVar.CurrentRunStepNumber=CtrlVar.CurrentRunStepNumber+1;
    
    if CtrlVar.PlotWaitBar 
        multiWaitbar('Run steps','Value',(CtrlVar.CurrentRunStepNumber-1-CtrlVar.CurrentRunStepNumber0)/CtrlVar.TotalNumberOfForwardRunSteps);
        multiWaitbar('Time','Value',CtrlVar.time/CtrlVar.TotalTime);
    end
    
    MUA=UpdateMUA(CtrlVar,MUA);
    
    % -adapt time step
    if CtrlVar.TimeDependentRun
        [CtrlVar.dt,dtRatio]=AdaptiveTimeStepping(CtrlVar,CtrlVar.time,CtrlVar.dt,RunInfo,dubdt,dvbdt,dhdt);
    end
    
    
    if CtrlVar.doDiagnostic
        if CtrlVar.InDiagnosticRunsDefineIceGeometryAtEveryRunStep
            [UserVar,s,b,S,B,alpha]=GetGeometry(UserVar,CtrlVar,MUA,CtrlVar.time,'sbSB');
            h=s-b;
        end
    elseif CtrlVar.DefineOceanSurfaceAtEachTimeStep
        [UserVar,~,~,S,~,~]=GetGeometry(UserVar,CtrlVar,MUA,CtrlVar.time,'S');
    end
    
    [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
    GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
    [UserVar,C,m]=GetSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
    [UserVar,AGlen,n]=GetAGlenDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
    if CtrlVar.UpdateBoundaryConditionsAtEachTimeStep
        [UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,CtrlVar.time,s,b,h,S,B,ub,vb,ud,vd,GF);
        [ub,vb,ud,vd]=StartVelocity(CtrlVar,MUA,BCs,ub,vb,ud,vd,s,b,h,S,B,rho,rhow,GF,AGlen,n,C,m);
    end
    
    
    %% adaptive meshing, adapt mesh, adapt-mesh
    if CtrlVar.AdaptMesh || CtrlVar.TimeGeometries.Flag ||  CtrlVar.FEmeshAdvanceRetreat
        
        % looks like far too many parameters here, but I need the whole model definition
        % and on return everything will be defined on a new mesh
        [UserVar,CtrlVar,MUA,BCs,GF,GLdescriptors,...
            s,b,h,S,B,ub,vb,ud,vd,l,rho,rhow,g,AGlen,n,C,m,ab,as,dhdt,dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1]=...
            AdaptMesh(UserVar,CtrlVar,MUA,BCs,...
            GF,GLdescriptors,alpha,...
            s,b,h,S,B,ub,vb,ud,vd,Ruv,Lubvb,l,rho,rhow,g,AGlen,n,C,m,ab,as,dhdt,dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1);
        
        
        if MUA.Nele==0 
            fprintf('FE mesh is empty \n ')
            break ;
        end
        
        if CtrlVar.AdaptMeshAndThenStop
            
            
            if CtrlVar.doplots  && CtrlVar.PlotMesh
                figure ; PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar)
            end
            
            save(CtrlVar.SaveInitialMeshFileName,'MUA') ;
            fprintf(CtrlVar.fidlog,' MUA was saved in %s .\n',CtrlVar.SaveInitialMeshFileName);
            fprintf('Exiting after remeshing because CtrlVar.AdaptMeshAndThenStop set to true. \n ')
            return
        end
        
    end
    
    %%
    
    ub0=ub ; vb0=vb; ud0=ud ; vd0=vd ; h0=h; s0=s ; b0=b;
    [UserVar,as0,ab0,dasdh0,dabdh0]=GetMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
    a0=as0+ab0;
    
    
    if ~CtrlVar.TimeDependentRun % Time independent run.  Solving for velocities for a given geometry (diagnostic steo).
        
        %% Diagnostic calculation (uv)
        if CtrlVar.InfoLevel >= 1 ; fprintf(CtrlVar.fidlog,' ==> Time independent step. Current run step: %i \n',CtrlVar.CurrentRunStepNumber) ;  end

        [UserVar,ub,vb,ud,vd,l,Kuv,Ruv,RunInfo,Lubvb]= uv(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF);

        dubdtm1=dubdt ; dvbdtm1=dvbdt; duddtm1=duddt ; dvddtm1=dvddt;
        
        if CtrlVar.dt==0 
            dubdt=zeros(MUA.Nnodes,1) ;  dvbdt=zeros(MUA.Nnodes,1);
            duddt=zeros(MUA.Nnodes,1) ;  dvddt=zeros(MUA.Nnodes,1);
        else
            dubdt=(ub-ub0)/CtrlVar.dt ; dvbdt=(vb-vb0)/CtrlVar.dt;
            duddt=(ud-ud0)/CtrlVar.dt ; dvddt=(vd-vd0)/CtrlVar.dt;
        end
        
        
    else   % Time-dependent run
        
        if CtrlVar.Implicituvh % Fully implicit time-dependent step (uvh)
            
            fprintf(CtrlVar.fidlog,...
                '\n ===== Implicit uvh going from t=%-.10g to t=%-.10g with dt=%-g. Done %-g %% of total time, and  %-g %% of steps. \n ',...
                CtrlVar.time,CtrlVar.time+CtrlVar.dt,CtrlVar.dt,100*CtrlVar.time/CtrlVar.TotalTime,100*(CtrlVar.CurrentRunStepNumber-1-CtrlVar.CurrentRunStepNumber0)/CtrlVar.TotalNumberOfForwardRunSteps);
            
            if CtrlVar.InitialDiagnosticStep   % if not a restart step, and if not explicitly requested by user, then do not do an inital dignostic step
                %% diagnostic step, solving for uv.  Always needed at a start of a transient run. Also done if asked by the user.
                RunInfo.converged=0; dIt=0;
                while RunInfo.converged==0
                    dIt=dIt+1;
                    
                    fprintf(CtrlVar.fidlog,' initial diagnostic step at t=%-.15g \n ',CtrlVar.time);
                    
                    [UserVar,ub,vb,ud,vd,l,Kuv,Ruv,RunInfo,Lubvb]= uv(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF);
                   
                    if  RunInfo.converged==1 || dIt==2
                        break
                    end
                    fprintf(CtrlVar.fidlog,' initial diagnostic step at t=%-g did not converge. Resetting (u,v) to zero and trying again \n ',CtrlVar.time);
                    TempFile=[CtrlVar.Experiment,'-UaDumpDiagnosticStepNotConverged.mat']; fprintf(CtrlVar.fidlog,' saving variables in %s \n ',TempFile) ; save(TempFile,'-v7.3')
                    ub=ub*0 ; vb=vb*0 ; ud=ud*0 ; vd=vd*0;  l.ubvb=l.ubvb*0; % if did not converge try to reset
                end
                
                if dIt==2 && RunInfo.converged==0
                    fprintf(CtrlVar.fidlog,' initial diagnostic step at t=%-g did not converge despite resetting (u,v) to zero \n ',CtrlVar.time);
                    fprintf(CtrlVar.fidlog,' saving variables in UaDump2 \n ') ; save UaDump2
                end
                
                ub0=ub ; ud0=ud ; vb0=vb ; vd0=vd;
                CtrlVar.InitialDiagnosticStep=0;
                
                if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 )
                    CtrlVar.UaOutputsInfostring='Diagnostic step';
                    CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
                    fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
                    UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,s,b,S,B,h,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,dasdh,dabdh,GF,BCs,l);
                 
                    if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls
                        fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
                            CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
                        return
                    end
                end
            end
            
            
            
            
            %% get an explicit estimate for u, v and h at the end of the time step
            
            [ub1,vb1,ud1,vd1,h1]=ExplicitEstimation(CtrlVar.dt,dtRatio,CtrlVar.CurrentRunStepNumber,ub,dubdt,dubdtm1,vb,dvbdt,dvbdtm1,ud,duddt,duddtm1,vd,dvddt,dvddtm1,h,dhdt,dhdtm1);
            
            %% advance the solution by dt using a fully implicit method with respect to u,v and h
            uvhStep=1;
            while uvhStep==1  && CtrlVar.dt > CtrlVar.dtmin % if uvh step does not converge, it is repeated with a smaller dt value
                
                [UserVar,as,ab,dasdh,dabdh]=GetMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time+CtrlVar.dt,s,b,h,S,B,rho,rhow,GF);
                as1=as ; ab1=ab;
                %        0  : values at t
                %        1  : explicit guess for values at t+dt
                % on output : converged values at t+dt
                
                
                [UserVar,ub,vb,ud,vd,h,l.ubvb,l.h,RunInfo,CtrlVar,BCs,CtrlVar.dt]=...
                    FIuvh2D(UserVar,CtrlVar,MUA,BCs,CtrlVar.dt,S,B,ub0,vb0,ud0,vd0,h0,ub1,vb1,ud1,vd1,h1,as0,ab0,as1,ab1,...
                    dubdt,dvbdt,duddt,dvddt,...
                    l.ubvb,l.h,...
                    AGlen,C,n,m,alpha,rho,rhow,g);
                uvhStep=0;  % assuming it converged
                
                if RunInfo.converged==0
                    
                    fprintf(CtrlVar.fidlog,' Warning : Reducing time step from %-g to %-g \n',CtrlVar.dt,CtrlVar.dt/10);
                    CtrlVar.dt=CtrlVar.dt/10;
                    uvhStep=1;  % continue within while loop
                    
                    fprintf(CtrlVar.fidlog,'Also resetting u1, v1, h1 to ub0, vb0 and h0, and setting estimates for Lagrange parameters to zero. \n');
                    [b,s,h0]=Calc_bs_From_hBS(h0,S,B,rho,rhow,CtrlVar,MUA.coordinates);
                    ub1=ub0*0 ; vb1=vb0*0 ;  ud1=ud0*0 ; vd1=vd0*0 ; h1=h0;
                    l.ubvb=l.ubvb*0; l.h=l.h*0;
                end
            end
            
            CtrlVar.time=CtrlVar.time+CtrlVar.dt;
            %CtrlVar.time=round(CtrlVar.time,14,'significant');
            
            [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
            
            dhdtm1=dhdt ; dubdtm1=dubdt ; dvbdtm1=dvbdt;
            if CtrlVar.dt==0 
                dhdt=zeros(MUA.Nnodes,1) ;  dubdt=zeros(MUA.Nnodes,1); dvbdt=zeros(MUA.Nnodes,1); dsdt=zeros(MUA.Nnodes,1) ; dbdt=zeros(MUA.Nnodes,1);
            else
                dhdt=(h-h0)/CtrlVar.dt; dubdt=(ub-ub0)/CtrlVar.dt ; dvbdt=(vb-vb0)/CtrlVar.dt; duddt=(ud-ud0)/CtrlVar.dt ; dvddt=(vd-vd0)/CtrlVar.dt; dsdt=(s-s0)/CtrlVar.dt; dbdt=(b-b0)/CtrlVar.dt;
            end
            
            
            
            
        elseif ~CtrlVar.Implicituvh % Semi-implicit time-dependent step. Implicit with respect to h, explicit with respect to u and v.
            
            if CtrlVar.InfoLevel>0 ; fprintf(CtrlVar.fidlog,'Semi-implicit transient step. Advancing time from t=%-g to t=%-g \n',CtrlVar.time,CtrlVar.time+CtrlVar.dt);end
            
            %% Diagnostic calculation (uv)
            if CtrlVar.InfoLevel >= 1 ; fprintf(CtrlVar.fidlog,' ==> Diagnostic step (uv). Current run step: %i \n',CtrlVar.CurrentRunStepNumber) ;  end
            tdiagnostic=tic;                  % -uv
            [UserVar,ub,vb,ud,vd,l,Kuv,Ruv,RunInfo,Lubvb]= uv(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF);
            
            
            tdiagnostic=toc(tdiagnostic);
            
            dubdtm1=dubdt ; dvbdtm1=dvbdt; duddtm1=duddt ; dvddtm1=dvddt;
            
            if CtrlVar.dt==0 
                dubdt=zeros(MUA.Nnodes,1) ;  dvbdt=zeros(MUA.Nnodes,1);
                duddt=zeros(MUA.Nnodes,1) ;  dvddt=zeros(MUA.Nnodes,1);
            else
                dubdt=(ub-ub0)/CtrlVar.dt ; dvbdt=(vb-vb0)/CtrlVar.dt;
                duddt=(ud-ud0)/CtrlVar.dt ; dvddt=(vd-vd0)/CtrlVar.dt;
            end
            
            tprognostic=tic;
            
            
            [ub1,vb1]=ExplicitEstimation(CtrlVar.dt,dtRatio,CtrlVar.CurrentRunStepNumber,ub,dubdt,dubdtm1,vb,dvbdt,dvbdtm1);
            
            dub1dt=dubdt; dvb1dt=dvbdt ;  dub0dt=dubdt; dvb0dt=dvbdt ; % could possibly be done a bit better
            [UserVar,as,ab,dasdh,dabdh]=GetMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time+CtrlVar.dt,s,b,h,S,B,rho,rhow,GF);
            a1=as+ab; da0dt=(a1-a0)/CtrlVar.dt ; da1dt=da0dt;
            if CtrlVar.dt==0 ; da0dt=zeros(MUA.Nnodes,1); da1dt=zeros(MUA.Nnodes,1) ; end
            [h,l]=SSS2dPrognostic(CtrlVar,MUA,BCs,l,h0,ub0,vb0,dub0dt,dvb0dt,a0,da0dt,ub1,vb1,a1,da1dt,dub1dt,dvb1dt);
            
            
            CtrlVar.time=CtrlVar.time+CtrlVar.dt;
            CtrlVar.time=round(CtrlVar.time,14,'significant');
            
            [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
            
            dhdtm1=dhdt ;
            if CtrlVar.dt==0 
                dhdt=zeros(MUA.Nnodes,1) ; dsdt=zeros(MUA.Nnodes,1) ; dbdt=zeros(MUA.Nnodes,1);
            else
                dhdt=(h-h0)/CtrlVar.dt; dsdt=(s-s0)/CtrlVar.dt; dbdt=(b-b0)/CtrlVar.dt;
            end
            
            
            %CtrlVar.hChange=1; CtrlVar.s=1; % h and s just changed
            tprognostic=toc(tprognostic);
            if CtrlVar.InfoLevel >= 1 && fprintf(CtrlVar.fidlog,'SSTREAM semi-implicit step in %-g sec, \t prognostic in %-g and diagnostic in %-g sec \n ',tprognostic+tdiagnostic,tprognostic,tdiagnostic) ; end
        end
    end
    
    %% calculations for this step are now done, only some plotting/writing issues do deal with
    
    % calculating derived quantities
    
    GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
    
    
    if CtrlVar.doPrognostic==1 && CtrlVar.InfoLevel>=10
        dhdtMean=(dhdt+dhdtm1)/2;
        [maxdhdt,imaxdhdt]=max(dhdtMean);
        [mindhdt,imindhdt]=min(dhdtMean);
        fprintf(CtrlVar.fidlog,'max(h) %-g \t min(h) %-g \t max(dhdt) %-g \t min(dhdt) %-g \t mean(dhdt) %-g \t median(dhdt) %-g \t rms(h) %-g \t h(max(dhdt)) %-g h(min(dhdt)) %-g\n ',...
            max(h),min(h),max(dhdtMean),min(dhdtMean),mean(dhdt),median(dhdt),sqrt(norm(dhdt)/numel(h)),h(imaxdhdt),h(imindhdt));
        
    end
    %% plotting results
    
    % UaOutputs
    
    if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 )
        CtrlVar.UaOutputsInfostring='inside transient loop and inside run-step loop';
        CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
        
        if CtrlVar.MassBalanceGeometryFeedback>0
             [UserVar,as,ab,dasdh,dabdh]=GetMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time+CtrlVar.dt,s,b,h,S,B,rho,rhow,GF);
        end
        
        fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
        
        UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,s,b,S,B,h,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,dasdh,dabdh,GF,BCs,l);
        if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls
            fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
                CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
            return
        end
    end
    
    if CtrlVar.WriteRestartFile==1 && mod(CtrlVar.CurrentRunStepNumber,CtrlVar.WriteRestartFileInterval)==0
        WriteRestartFile()
    end
    
    %% write a dump file
    
    if CtrlVar.WriteDumpFile && (ReminderFraction(CtrlVar.time,CtrlVar.WriteDumpFileTimeInterval)<1e-5 || ReminderFraction(CtrlVar.CurrentRunStepNumber,CtrlVar.WriteDumpFileStepInterval)==0)
        
        
        if ReminderFraction(CtrlVar.time,CtrlVar.WriteDumpFileTimeInterval)<1e-5
            dumpfile=sprintf('%s-DumpFile%-gT-',CtrlVar.UserVar,CtrlVar.time);  dumpfile=regexprep(dumpfile,'\.','k');
        else
            dumpfile=sprintf('%sDumpFile%i',CtrlVar.UserVar,CtrlVar.CurrentRunStepNumber);
        end
        
        fprintf(CtrlVar.fidlog,' Saving everything in file %s  at t=%-g \n',dumpfile,CtrlVar.time);
        
        try
            save(dumpfile','s','h','b','B','S','ub','vb','wSurf','wBed','as','ab','time','MUA','AGlen','C',...
                'Luv','Luvrhs','lambdauv','Lh','Lhrhs','lambdah','n','m','alpha','rhow','rho','g',...
                'dubdt','dvbdt','a0','da0dt','dhdt','dhdtm1','dubdtm1','dvbdtm1','DTxy','TRIxy','CtrlVar','-v7.3')
        catch exception
            fprintf(CtrlVar.fidlog,' Could not save file %s \n ',dumpfile);
            fprintf(CtrlVar.fidlog,'%s \n',exception.message);
        end
    end
    
end

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
        [UserVar,as,ab,dasdh,dabdh]=GetMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time+CtrlVar.dt,s,b,h,S,B,rho,rhow,GF);
    end
    
    fprintf(' Calling UaOutputs. UaOutputsInfostring=%s , UaOutputsCounter=%i \n ',CtrlVar.UaOutputsInfostring,CtrlVar.UaOutputsCounter)
    
    UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,s,b,S,B,h,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,dasdh,dabdh,GF,BCs,l);
    if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls
        fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
            CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
        return
    end
end



%% saving outputs

if CtrlVar.WriteRestartFile==1 &&  mod(CtrlVar.CurrentRunStepNumber,CtrlVar.WriteRestartFileInterval)~=0
    
    WriteRestartFile()
    
end

if CtrlVar.PlotWaitBar ;     multiWaitbar('CloseAll'); end
tTime=toc(tTime); fprintf(CtrlVar.fidlog,' Total time : %-g sec \n',tTime) ;


if CtrlVar.fidlog~= 1 ; fclose(CtrlVar.fidlog); end
fclose(CtrlVar.InfoFile);

SayGoodbye(CtrlVar)

%% nested functions


    function WriteRestartFile
        RestartFile=CtrlVar.NameOfRestartFiletoWrite;
        fprintf(CtrlVar.fidlog,' \n ################## %s %s ################### \n Writing restart file %s  at t=%-g \n %s \n ',CtrlVar.Experiment,datestr(now),RestartFile,CtrlVar.time);
        %[DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);
        CtrlVarInRestartFile=CtrlVar;
        nStep=CtrlVar.CurrentRunStepNumber;  % later get rid of nStep from all restart files
        Itime=CtrlVar.CurrentRunStepNumber;  % later get rid of Itime from all restart files
        time=CtrlVar.time;
        dt=CtrlVar.dt;
        try
            save(RestartFile,'CtrlVarInRestartFile','UserVar','MUA','BCs','time','dt','s','b','S','B','h','ub','vb','ud','vd','dhdt','dsdt','dbdt','C','AGlen','m','n','rho','rhow','as','ab','GF',...
                'nStep','Itime','dhdtm1','dubdt','dvbdt','dubdtm1','dvbdtm1','duddt','dvddt','duddtm1','dvddtm1',...
                'GLdescriptors','l','alpha','g','-v7.3');
            
            
            fprintf(CtrlVar.fidlog,' Writing restart file was successful. \n');
            
        catch exception
            fprintf(CtrlVar.fidlog,' Could not save restart file %s \n ',RestartFile);
            fprintf(CtrlVar.fidlog,'%s \n',exception.message);
        end
        
    end



    function WriteAdjointRestartFile
        
        
        fprintf(CtrlVar.fidlog,'Saving adjoint restart file: %s \n ',CtrlVar.NameOfAdjointRestartFiletoWrite);
        save(CtrlVar.NameOfAdjointRestartFiletoWrite,...
            'UserVar','CtrlVar','MUA','BCs','s','b','h','S','B','ub','vb','ud','vd','alpha','rho','rhow','g',...
            'as','ab','GF','BCs','l',...
            'InvStartValues','Priors','Meas','BCsAdjoint','Info','InvFinalValues','xAdjoint','yAdjoint','-v7.3');
        
        
        if CtrlVar.AGlenisElementBased
            xA=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1));
            yA=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2));
        else
            xA=MUA.coordinates(:,1);
            yA=MUA.coordinates(:,2);
        end
        
        if CtrlVar.CisElementBased
            xC=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1));
            yC=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2));
        else
            xC=MUA.coordinates(:,1);
            yC=MUA.coordinates(:,2);
        end
        
        fprintf(CtrlVar.fidlog,' saving C and m  in file %s \n ',CtrlVar.NameOfFileForSavingSlipperinessEstimate)        ;
        save(CtrlVar.NameOfFileForSavingSlipperinessEstimate,'C','m','xC','yC','MUA')
        
        fprintf(CtrlVar.fidlog,' saving AGlen and m in file %s \n ',CtrlVar.NameOfFileForSavingAGlenEstimate) ;
        save(CtrlVar.NameOfFileForSavingAGlenEstimate,'AGlen','n','xA','yA','MUA')
        
    end
































end
