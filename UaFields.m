classdef UaFields
    
  
    properties
        
        x=[];
        y=[];
        
        ub=[];
        vb=[];
        
        ud=[];
        vd=[];
        
        uo=[];
        vo=[];
        
        ua=[];
        va=[];
        
        
        
        s=[] ;
        sInit=[]; 
        
        b=[];
        bmin=[];
        bmax=[];
        bInit=[]
        
        h=[];
        hInit=[]; 
        S=[];
        
        B=[];
        Bmin=[];
        Bmax=[];
        BInit=[];
        
        AGlen=[];
        AGlenmin=[];
        AGlenmax=[];
        
        C=[];
        Cmin=[];
        Cmax=[];
        m=[];
        n=[];
        rho=[];
        rhow=[];
        
        q=[]; 
        muk=[];
        
        Co=[];
        mo=[]
        Ca=[];
        ma=[];
        
        as=[];
        ab=[];
        dasdh=[];
        dabdh=[];

        
        dhdt=[] ;
        dsdt=[] ;
        dbdt=[] ;
     
        dubdt=[];
        dvbdt=[];
        
        duddt=[];
        dvddt=[];

      
        
        g=[];
        alpha=[];
        
        time=[]; 
        dt=[] ; 
        
        GF=[];
        GFInit=[];
        
        LSF=[] % Level Set Field
        LSFMask=[]; 
        LSFnodes=[];
        c=[] ; % calving rate
        LSFqx=[] ;
        LSFqy=[] ; 

        N=[] ; 
        aw=[];
        hw=[]; 
        phi=[] ;
        uw=[];
        vw=[];

        txx=[];
        txy=[];
        tyy=[];

      
    end
    
end