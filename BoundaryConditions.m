classdef BoundaryConditions
    %%
    %
    % Very simple class specifying the different types of Dirichlet boundary conditions.
    %
    %
    %
    %%
    
    properties
        
        % basal velocities
        ubFixedNode=[];
        ubFixedValue=[];
        vbFixedNode=[];
        vbFixedValue=[];
        
        ubTiedNodeA=[];
        ubTiedNodeB=[];
        vbTiedNodeA=[];
        vbTiedNodeB=[];
        
        ubvbFixedNormalNode=[];
        ubvbFixedNormalValue=[];
        
        % deformational velocity, only used when solving SIA or hybrid
        udFixedNode=[];
        udFixedValue=[];
        vdFixedNode=[];
        vdFixedValue=[];
        
        udTiedNodeA=[];
        udTiedNodeB=[];
        vdTiedNodeA=[];
        vdTiedNodeB=[];
        
        udvdFixedNormalNode=[];
        udvdFixedNormalValue=[];
        
        % thickness  - only used if solving for thickness, for example in a uvh step
        hFixedNode=[];
        hFixedValue=[];
        hTiedNodeA=[];
        hTiedNodeB=[];
        
        % positive thickness constraints - these are introduced automatically internally, do not
        % prescribe directly! 
        hPosNode=[];
        hPosValue=[];
        
        % Boundary conditions for the level set field
        % Only need when using the level-set method (currenlty under development, so don't
        % use!)
        LSFFixedNode=[];
        LSFFixedValue=[];
        LSFTiedNodeA=[];
        LSFTiedNodeB=[];
        LSFL=[] ; 
        LSFrhs=[] ; 

        % rate of thickness change - only used when using calculating dh/dt in combination with the ajoint methods
        %                            for example when using measurements of dh/dt in an inversion
        dhdtFixedNode=[];
        dhdtFixedValue=[];
        dhdtTiedNodeA=[];
        dhdtTiedNodeB=[];
        
        
    end
    
end