classdef UaRunInfo
    
    properties
        
        Inverse
        Forward 
    end
    
    
    methods
        
        function obj = UaRunInfo()
            
            
            obj.Inverse.J = NaN ;
            obj.Inverse.Iterations = 0;
            obj.Inverse.J = NaN ;
            obj.Inverse.I = NaN ;
            obj.Inverse.R = NaN ;
            obj.Inverse.StepSize = NaN ;
            
            
            obj.Forward.Converged=0;
            obj.Forward.Iterations=NaN;
            obj.Forward.Residual=NaN;
            
            
        end
    end
    
end
