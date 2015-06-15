function [rho,rhow,g]=DefineDensities(Experiment,CtrlVar,MUA,time,s,b,h,S,B)

    
    warning('Ua:DefaultDefine','Using default DefineDensities \n')
        
    rho=900+zeros(MUA.Nnodes,1) ; rhow=1030; g=9.81/1000;
    
    
end
