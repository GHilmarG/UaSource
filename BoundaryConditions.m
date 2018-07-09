classdef BoundaryConditions
    %%
    % 
    %
    %
    %%
    
    properties
        
        % basal
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
        
        % deformational
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
        
        % thickness
        hFixedNode=[];
        hFixedValue=[];
        hTiedNodeA=[];
        hTiedNodeB=[];
        
        hPosNode=[];
        hPosValue=[];
         
    end
 
end