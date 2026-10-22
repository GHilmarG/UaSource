classdef UaLagrangeVariables
    
    properties
        
        ubvb=[];
        udvd=[];
        h=[];
        dth=NaN;   % (9 Oct 2026) time step for which the multipliers of the thickness constraints (h) were calculated.
                   % These multipliers scale with the time step and are rescaled when used for a different time step.
    end
    
end