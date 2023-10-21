classdef PriorProbabilityDistribution
    
    properties
        
        
        B=[];
        Bmax=[];
        Bmin=[];
        CovB=[];
        
        h=[]; 
        
        rho=[];
        rhow=[];
        m=[];
        n=[];
        q=[];
        muk=[];
        
        AGlen=[];
        CovAGlen=[]
        AGlenmax=[];
        AGlenmin=[];
        
        C=[];
        CovC=[];
        Cmax=[];
        Cmin=[];
        
        % These fields are for saving the `true' values, generally only makes sense if doing some synthetic inversion
        % experiments, in which case these fields are those used to generate the synthetic measurements. 
        TrueC=[];
        TrueAGlen=[];
        Trueb=[];
        TrueB=[];
        
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
            
            obj.Regularize.B.gs=[];
            obj.Regularize.B.ga=[];
        end
    end
    
end