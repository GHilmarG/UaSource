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
        
        %         Cgs=1;
        %         Cga=1;
        %         logCga=1;
        %         logCgs=1 ;
        %
        %         AGlengs=1;
        %         AGlenga=1;
        %         logAGlenga=1;
        %         logAGlengs=1 ;
        %
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