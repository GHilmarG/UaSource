classdef PriorProbabilityDistribution
    
    properties
        
        s=[];
        b=[];
        S=[];
        B=[];
        rho=[];
        rhow=[];
        m=[];
        n=[];
        
        C=[];
        AGlen=[];
        CovAGlen=[]
        CovC=[];
        
        TrueC=[];
        TrueAGlen=[];
 
        Regularize
        
        
    end
    
    methods
        
        function obj=PriorProbabilityDistribution()
            
            obj.Regularize.C.gs=[];
            obj.Regularize.C.ga=[];
            obj.Regularize.logC.gs=[];
            obj.Regularize.logC.ga=[];
            
            obj.Regularize.AGlen.gs=[];
            obj.Regularize.AGlen.ga=[];
            obj.Regularize.logAGlen.gs=[];
            obj.Regularize.logAGlen.ga=[];
            
            
        end
    end
    
end