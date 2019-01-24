classdef PriorProbabilityDistribution
    
    properties
        
       

        b=[];
        Covb=[];
        bmax=[];
        bmin=[];
        
        B=[];
        Bmax=[];
        Bmin=[];
        CovB=[];
        
        
        rho=[];
        rhow=[];
        m=[];
        n=[];
        
        AGlen=[];
        CovAGlen=[]
        AGlenmax=[];
        AGlenmin=[];
        
        C=[];
        CovC=[];
        Cmax=[];
        Cmin=[];
        
        TrueC=[];
        TrueAGlen=[];
        Trueb=[];
        
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
            
            obj.Regularize.b.gs=[];
            obj.Regularize.b.ga=[];
            
            obj.Regularize.B.gs=[];
            obj.Regularize.B.ga=[];
        end
    end
    
end