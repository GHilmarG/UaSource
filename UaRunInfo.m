classdef UaRunInfo
    
    properties
        
        Inverse
        Forward
        BackTrack
        CPU
    end
    
    
    methods
        
        function obj = UaRunInfo()
            
            
            obj.Inverse.J = NaN ;
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
            obj.Forward.Residual=NaN;
            
            obj.Forward.ActiveSetConverged=NaN;
            
            obj.BackTrack.Converged=NaN;
            obj.BackTrack.iarm=NaN;
            obj.BackTrack.Infovector=NaN;
            
            obj.CPU.Assembly=0;
            obj.CPU.Solution=0;
            
            
        end
    end
    
end
