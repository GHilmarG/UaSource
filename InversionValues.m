
classdef InversionValues
    
    properties
        
        
        m=[];
        n=[];
        q=[];
        muk=[];
        C=[];
        AGlen=[];
        B=[];
        h=[] ;  % Only use this for storing variables 
        
        SearchStepSize=[];
        uAdjoint=[];
        vAdjoint=[];
        
        J=[];
        I=[];
        R=[];
        
        RAGlen=[];
        RC=[];

        RAa=[] ;
        RAs=[] ;
        RCa=[] ;
        RCs=[] ;
        
        dJdp=[];
        dJdC=[];
        dJdAGlen=[];
        dJdb=[]; 
        dJdB=[]; 
        
        dIdp=[];
        dIdC=[];
        dIdAGlen=[];
        dIdb=[]; 
        
        dRdp=[];
        dRdC=[];
        dRdAGlen=[];
        dRb=[]; 
        
        dJdpTest=[];
        dJdCTest=[];
        dJdAGlenTest=[];
        dJdbTest=[]; 
        dJdBTest=[]; 
        
        dIdpTest=[];
        dIdCTest=[];
        dIdAGlenTest=[];
        dIdbTest=[]; 
        dIdBTest=[]; 
    
        
    end
    
end