classdef UaRunInfo
    
    properties
        
        Inverse
        Forward
        BackTrack
        CPU
        Message
        MeshAdapt
        File
        
        
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
            obj.Inverse.GradNorm = NaN ;
            obj.Inverse.ConjGradUpdate = 0 ;
            obj.Inverse.fmincon=struct;
            obj.Inverse.fminunc=struct;
            
            obj.Forward.Converged=0;
            obj.Forward.Iterations=NaN;
            obj.Forward.IterationsTotal=0;
            obj.Forward.Residual=NaN;
            obj.Forward.time=NaN;
            obj.Forward.dt=NaN;
            
            
            obj.Forward.ActiveSetConverged=NaN;
            
            obj.BackTrack.Converged=NaN;
            obj.BackTrack.iarm=NaN;
            obj.BackTrack.Infovector=NaN;
            
            
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
            
            
        end
        
    end
    
end
