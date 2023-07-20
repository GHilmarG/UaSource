classdef (ConstructOnLoad) UaRunInfo

    properties

        Inverse
        Forward
        BackTrack
        CPU
        Message
        MeshAdapt
        File
        Mapping
        LevelSet
    end


    methods

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

            N=1000; % initial memory allocation
            obj.Forward.uvConverged=false;
            obj.Forward.uvIterations=NaN(N,1);
            obj.Forward.uvResidual=NaN(N,1);
            obj.Forward.uvBackTrackSteps=NaN(N,1);



            obj.Forward.dtRestart=NaN;


            N=1000; % initial memory allocation
            obj.Forward.time=zeros(N,1)+NaN;
            obj.Forward.dt=zeros(N,1)+NaN;
            obj.Forward.uvhConverged=false;
            obj.Forward.uvhIterations=zeros(N,1)+NaN;
            obj.Forward.uvhResidual=zeros(N,1)+NaN;
            obj.Forward.uvhBackTrackSteps=zeros(N,1)+NaN;
            obj.Forward.uvhActiveSetIterations=zeros(N,1)+NaN;
            obj.Forward.uvhActiveSetCyclical=zeros(N,1)+NaN;
            obj.Forward.uvhActiveSetConstraints=zeros(N,1)+NaN;

            obj.Forward.hConverged=0;
            obj.Forward.hIterations=NaN(N,1);
            obj.Forward.hResidual=NaN(N,1);
            obj.Forward.hBackTrackSteps=zeros(N,1)+NaN;
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
            obj.LevelSet.time=zeros(N,1)+NaN;
            obj.LevelSet.Iterations=zeros(N,1)+NaN;
            obj.LevelSet.Residual=zeros(N,1)+NaN;
            obj.LevelSet.BackTrackSteps=zeros(N,1)+NaN;
            obj.LevelSet.Phase=strings(N,1);

            obj.CPU.Total=0;
            obj.CPU.Assembly.uv=0;
            obj.CPU.Solution.uv=0;
            obj.CPU.Assembly.uvh=0;
            obj.CPU.Solution.uvh=0;
            obj.CPU.WallTime="";


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

    end

    methods (Static)
        function obj = loadobj(s)

            obj=s;

            % Make sure the loaded UaRunInfo is up-to date
            % add in here any new modifications
            if ~isfield(s.Forward,'AdaptiveTimeSteppingResetCounter')
                obj.Forward.AdaptiveTimeSteppingResetCounter=0;
            end

            if ~isfield(s.Forward,'uvhIterations')

                N=1000; % initial memory allocation
                obj.Forward.time=zeros(N,1)+NaN;
                obj.Forward.dt=zeros(N,1)+NaN;
                obj.Forward.uvhIterations=zeros(N,1)+NaN;
                obj.Forward.uvhResidual=zeros(N,1)+NaN;
                obj.Forward.uvhBackTrackSteps=zeros(N,1)+NaN;
                obj.Forward.uvhActiveSetIterations=zeros(N,1)+NaN;
                obj.Forward.uvhActiveSetCyclical=zeros(N,1)+NaN;
                obj.Forward.uvhActiveSetConstraints=zeros(N,1)+NaN;

            end

            if ~isfield(s.BackTrack,'iarm')

                obj.BackTrack.Converged=NaN;
                obj.BackTrack.iarm=NaN;
                obj.BackTrack.Infovector=NaN;
                obj.BackTrack.nFuncEval=NaN;
                obj.BackTrack.nExtrapolationSteps=NaN;


            end

            if ~isfield(s.Forward,'hiCount')
                obj.Forward.hiCount=0;
            end

            %             I got rid of this counter, if I keep this in here, these fields will always be set to zero on load
            %             if ~isfield(s.Forward,'iCounter')
            %
            %                 N=1000; % initial memory allocation
            %                 obj.Forward.time=zeros(N,1)+NaN;
            %                 obj.Forward.dt=zeros(N,1)+NaN;
            %                 obj.Forward.uvhIterations=zeros(N,1)+NaN;
            %                 obj.Forward.uvhResidual=zeros(N,1)+NaN;
            %                 obj.Forward.uvhBackTrackSteps=zeros(N,1)+NaN;
            %                 obj.Forward.uvhActiveSetIterations=zeros(N,1)+NaN;
            %                 obj.Forward.uvhActiveSetCyclical=zeros(N,1)+NaN;
            %                 obj.Forward.uvhActiveSetConstraints=zeros(N,1)+NaN;
            %             end

            if ~isfield(s.Forward,'ubvbRecalculatedOnNewMesh')
                obj.Forward.ubvbRecalculatedOnNewMesh=false;
            end

            if ~isfield(s.LevelSet,'SolverConverged')
                obj.LevelSet.SolverConverged=0;
                obj.LevelSet.Iterations=NaN;
                obj.LevelSet.rForce=NaN;
                obj.LevelSet.rWork=NaN;
            end

            if ~isfield(s.Inverse,'nFuncEval')
                obj.Inverse.nFuncEval = 0 ;
            end


        end

    end

end
