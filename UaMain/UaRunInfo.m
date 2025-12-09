

classdef (ConstructOnLoad) UaRunInfo < matlab.mixin.CustomElementSerialization 

    properties

        Inverse
        Forward
        BackTrack
        CPU
        WallTime
        Message
        MeshAdapt
        File
        Mapping
        LevelSet
    end


    methods

        % Constructor
        function obj = UaRunInfo()


            obj.File.fid = NaN ;
            obj.File.Name = NaN ;


            obj.Inverse.Iterations = 0;
            obj.Inverse.J = NaN ;
            obj.Inverse.I = NaN ;
            obj.Inverse.R = NaN ;
            obj.Inverse.StepSize = NaN ;
            obj.Inverse.GradNorm = NaN ;
            obj.Inverse.nFuncEval = 0 ;
            obj.Inverse.ConjGradUpdate = 0 ;
            obj.Inverse.fmincon=struct;
            obj.Inverse.fminunc=struct;

            N=2; % initial memory allocation
            obj.Forward.uvConverged=false;
            obj.Forward.uvIterations=NaN(N,1);
            obj.Forward.uvResidual=NaN(N,1);
            obj.Forward.uvBackTrackSteps=NaN(N,1);

            obj.Forward.uv2hIterations=NaN(N,1);

            obj.Forward.dtRestart=NaN;


            N=2; % initial memory allocation
            obj.Forward.time=NaN(N,1);
            obj.Forward.dt=NaN(N,1);
            obj.Forward.uvhConverged=false;
            obj.Forward.uvhIterations=NaN(N,1);
            obj.Forward.uvhResidual=NaN(N,1);
            obj.Forward.uvhBackTrackSteps=NaN(N,1);

            obj.Forward.uvhActiveSetIterations=NaN(N,1);
            obj.Forward.uvhActiveSetCyclical=NaN(N,1);
            obj.Forward.uvhActiveSetConstraints=NaN(N,1);

            obj.Forward.hActiveSetIterations=NaN(N,1);
            obj.Forward.hActiveSetCyclical=NaN(N,1);
            obj.Forward.hActiveSetConstraints=NaN(N,1);

            obj.Forward.hConverged=0;
            obj.Forward.hIterations=NaN(N,1);
            obj.Forward.hResidual=NaN(N,1);
            obj.Forward.hBackTrackSteps=NaN(N,1);
            obj.Forward.hiCount=0;


            obj.Forward.ubvbRecalculatedOnNewMesh=false;

            obj.Forward.ActiveSetConverged=NaN;

            obj.Forward.AdaptiveTimeSteppingResetCounter=0;
            obj.Forward.AdaptiveTimeSteppingTimeStepModifiedForOutputs=0;


            obj.BackTrack.Converged=NaN;
            obj.BackTrack.iarm=0;
            obj.BackTrack.Infovector=NaN;
            obj.BackTrack.nFuncEval=NaN;
            obj.BackTrack.nExtrapolationSteps=NaN;

            obj.LevelSet.iCount=0;
            obj.LevelSet.time=NaN(N,1);
            obj.LevelSet.Iterations=NaN(N,1);
            obj.LevelSet.Residual=NaN(N,1);
            obj.LevelSet.BackTrackSteps=NaN(N,1);
            obj.LevelSet.Phase=strings(N,1);


            obj.CPU.tic=duration(0,0,0);
            obj.CPU.toc=duration(0,0,0);
          
            obj.CPU.Total=duration(0,0,0);
            

            obj.CPU.Assembly.uv=[];
            obj.CPU.Solution.uv=[];
            obj.CPU.Assembly.uvh=[];
            obj.CPU.Solution.uvh=[];

            obj.CPU.Inversion=[];

            obj.WallTime.Total=duration(0,0,0);
            obj.WallTime.tic=duration(0,0,0);
            obj.WallTime.toc=duration(0,0,0);

            obj.Message="" ;

            obj.MeshAdapt.Method="";
            obj.MeshAdapt.isChanged=false;
            obj.MeshAdapt.Mesh.Nele=NaN;
            obj.MeshAdapt.Mesh.Nnodes=NaN;
            obj.MeshAdapt.Mesh.RunStepNumber=NaN;
            obj.MeshAdapt.Mesh.time=NaN;

            obj.Mapping.nNewNodes=NaN;
            obj.Mapping.nOldNodes=NaN;
            obj.Mapping.nIdenticalNodes=NaN;
            obj.Mapping.nNotIdenticalNodes=NaN;
            obj.Mapping.nNotIdenticalNodesOutside=NaN;
            obj.Mapping.nNotIdenticalNodesInside=NaN;


        end
        % 
        % function  obj=set.Forward(obj,val)
        % 
        %     obj.Forward=val ;
        % 
        % 
        % end
        % 
        % function  a=get.Forward(obj)
        % 
        %     fprintf("here\n")
        %     a=obj.Forward;
        % 
        % 
        % end


        function obj=ExtendAllocations(obj,CtrlVar,index)

            nElements=numel(obj.Forward.time) ;

            if index>nElements
                nPadding=max(2*nElements,index) ;
                Padding=nan(nPadding,1);
            else
                return
            end

            if contains(CtrlVar.ForwardTimeIntegration,["-uvh-","-uv-"])



                obj.Forward.time=UaResize(obj.Forward.time,nPadding,FillValue=nan) ; 
                obj.Forward.dt=UaResize(obj.Forward.dt,nPadding,FillValue=nan) ; 

                obj.Forward.uvhIterations=UaResize(obj.Forward.uvhIterations,nPadding,FillValue=nan) ; 
                obj.Forward.uvhResidual=UaResize(obj.Forward.uvhResidual,nPadding,FillValue=nan) ; 
                obj.Forward.uvhBackTrackSteps=UaResize(obj.Forward.uvhBackTrackSteps,nPadding,FillValue=nan) ; 

                obj.Forward.uvhActiveSetIterations=UaResize(obj.Forward.uvhActiveSetIterations,nPadding,FillValue=nan) ; 
                obj.Forward.uvhActiveSetCyclical=UaResize(obj.Forward.uvhActiveSetCyclical,nPadding,FillValue=nan) ; 
                obj.Forward.uvhActiveSetConstraints=UaResize(obj.Forward.uvhActiveSetConstraints,nPadding,FillValue=nan) ; 


                obj.Forward.hActiveSetIterations=UaResize(obj.Forward.hActiveSetIterations,nPadding,FillValue=nan) ; 
                obj.Forward.hActiveSetCyclical=UaResize(obj.Forward.hActiveSetCyclical,nPadding,FillValue=nan) ; 
                obj.Forward.hActiveSetConstraints=UaResize(obj.Forward.hActiveSetConstraints,nPadding,FillValue=nan) ; 

                obj.Forward.uvIterations=UaResize(obj.Forward.uvIterations,nPadding,FillValue=nan) ; 
                obj.Forward.uvResidual=UaResize(obj.Forward.uvResidual,nPadding,FillValue=nan) ; 
                obj.Forward.uvBackTrackSteps=UaResize(obj.Forward.uvBackTrackSteps,nPadding,FillValue=nan) ; 

                obj.Forward.uv2hIterations=UaResize(obj.Forward.uv2hIterations,nPadding,FillValue=nan) ; 


            end

            if contains(CtrlVar.ForwardTimeIntegration,"-h-")

                obj.Forward.hIterations=[obj.Forward.hIterations;Padding];
                obj.Forward.hResidual=[obj.Forward.hResidual;Padding];
                obj.Forward.hBackTrackSteps=[obj.Forward.hBackTrackSteps;Padding];


            end

            if CtrlVar.LevelSetMethod

                Padding=NaN(nPadding,1);
                obj.LevelSet.time=[obj.LevelSet.time;Padding];
                obj.LevelSet.Iterations=[obj.LevelSet.Iterations;Padding];
                obj.LevelSet.Residual=[obj.LevelSet.Residual;Padding];
                obj.LevelSet.BackTrackSteps=[obj.LevelSet.BackTrackSteps;Padding];
                obj.LevelSet.Phase=[obj.LevelSet.Phase;strings(nPadding,1)];


            end
        end
    end

    methods (Static)
        % function obj = loadobj(s)
        % 
        %     obj=s;
        % 
        %     % Make sure the loaded UaRunInfo is up-to date
        %     % add in here any new modifications
        %     if ~isfield(s.Forward,'AdaptiveTimeSteppingResetCounter')
        %         obj.Forward.AdaptiveTimeSteppingResetCounter=0;
        %     end
        % 
        %     if ~isfield(s.Forward,'uvhIterations')
        % 
        %         N=2; % initial memory allocation
        %         obj.Forward.time=NaN(N,1);
        %         obj.Forward.dt=NaN(N,1);
        %         obj.Forward.uvhIterations=NaN(N,1);
        %         obj.Forward.uvhResidual=NaN(N,1);
        %         obj.Forward.uvhBackTrackSteps=NaN(N,1);
        % 
        %         obj.Forward.uvhActiveSetIterations=NaN(N,1);
        %         obj.Forward.uvhActiveSetCyclical=NaN(N,1);
        %         obj.Forward.uvhActiveSetConstraints=NaN(N,1);
        % 
        %     end
        % 
        %     if ~isfield(s.Forward,'hActiveSetIterations')
        %         N=2; 
        %         obj.Forward.hActiveSetIterations=NaN(N,1);
        %         obj.Forward.hActiveSetCyclical=NaN(N,1);
        %         obj.Forward.hActiveSetConstraints=NaN(N,1);
        % 
        %     end
        % 
        %     if ~isfield(s.BackTrack,'iarm')
        % 
        %         obj.BackTrack.Converged=NaN;
        %         obj.BackTrack.iarm=NaN;
        %         obj.BackTrack.Infovector=NaN;
        %         obj.BackTrack.nFuncEval=NaN;
        %         obj.BackTrack.nExtrapolationSteps=NaN;
        % 
        % 
        %     end
        % 
        %     if ~isfield(s.Forward,'hiCount')
        %         obj.Forward.hiCount=0;
        %     end
        % 
        % 
        %     if ~isfield(s.Forward,'ubvbRecalculatedOnNewMesh')
        %         obj.Forward.ubvbRecalculatedOnNewMesh=false;
        %     end
        % 
        %     if ~isfield(s.LevelSet,'SolverConverged')
        %         obj.LevelSet.SolverConverged=0;
        %         obj.LevelSet.Iterations=NaN;
        %         obj.LevelSet.rForce=NaN;
        %         obj.LevelSet.rWork=NaN;
        %     end
        % 
        %     if ~isfield(s.Inverse,'nFuncEval')
        %         obj.Inverse.nFuncEval = 0 ;
        %     end
        %
        %
        %     if ~isfield(s.CPU,'tic')
        %         s.CPU.tic=0;
        %     end
        % end


        function modifyIncomingSerializationContent(sObj)


            % Make sure the loaded UaRunInfo is up-to date
            % add in here any new modifications
            if ~isfield(sObj.Forward,'AdaptiveTimeSteppingResetCounter')
                sObj.Forward.AdaptiveTimeSteppingResetCounter=0;
            end

            if ~isfield(sObj.Forward,'uvhIterations')

                N=2; % initial memory allocation
                sObj.Forward.time=NaN(N,1);
                sObj.Forward.dt=NaN(N,1);
                sObj.Forward.uvhIterations=NaN(N,1);
                sObj.Forward.uvhResidual=NaN(N,1);
                sObj.Forward.uvhBackTrackSteps=NaN(N,1);

                sObj.Forward.uvhActiveSetIterations=NaN(N,1);
                sObj.Forward.uvhActiveSetCyclical=NaN(N,1);
                sObj.Forward.uvhActiveSetConstraints=NaN(N,1);

            end

            if ~isfield(sObj.Forward,'hActiveSetIterations')
                N=2; 
                sObj.Forward.hActiveSetIterations=NaN(N,1);
                sObj.Forward.hActiveSetCyclical=NaN(N,1);
                sObj.Forward.hActiveSetConstraints=NaN(N,1);

            end

            if ~isfield(sObj.BackTrack,'iarm')

                sObj.BackTrack.Converged=NaN;
                sObj.BackTrack.iarm=NaN;
                sObj.BackTrack.Infovector=NaN;
                sObj.BackTrack.nFuncEval=NaN;
                sObj.BackTrack.nExtrapolationSteps=NaN;


            end

            if ~isfield(sObj.Forward,'hiCount')
                sObj.Forward.hiCount=0;
            end


            if ~isfield(sObj.Forward,'ubvbRecalculatedOnNewMesh')
                sObj.Forward.ubvbRecalculatedOnNewMesh=false;
            end

     

            if ~isfield(sObj.Inverse,'nFuncEval')
                sObj.Inverse.nFuncEval = 0 ;
            end





            if ~isfield(sObj.CPU,"tic")
                sObj.CPU.tic=duration(0,0,0);
                sObj.CPU.toc=duration(0,0,0);
            end

            if ~isfield(sObj.CPU,"Total")
                sObj.CPU.Total=duration(0,0,0);
            end

            if ~sObj.hasNameValue("WallTime")
               sObj.addNameValue("WallTime",struct("Total",duration(0,0,0),"tic",duration(0,0,0),"toc",duration(0,0,0)))
            end
     

            if ~isfield(sObj.WallTime,"Total")
                sObj.WallTime.Total=duration(0,0,0);
                sObj.WallTime.tic=duration(0,0,0);
                sObj.WallTime.toc=duration(0,0,0);
            end

            if ~sObj.hasNameValue("LevelSet")

                N=2;
                sObj.addNameValue("LevelSet",struct("iCount",0,"time",NaN(N,1),"Iterations",NaN(N,1),"Residual",NaN(N,1),"BackTrackSteps",NaN(N,1),"Phase",strings(N,1)))
           
            end

       

            %if sObj.hasNameValue("Name")
            %     nameArray = split(sObj.Name);
            %     sObj.addNameValue("FirstName",nameArray(1));
            %     sObj.addNameValue("LastName",nameArray(2));
            %     sObj.remove("Name");
            %     sObj.rename("Department","Division");
            %     id = split(sObj.EmployeeID,"-");
            %     sObj.updateValue("EmployeeID",id(2));
            % end
            %end

        end
    end

end
