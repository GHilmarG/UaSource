function [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1)
    
  
    narginchk(15,15)
    nargoutchk(2,9) 
    
    if ~(nargout==2 || nargout==9)
        error('af')
    end
    
    
    
    switch lower(CtrlVar.LevelSetAssembly)
        
        case "consistent"
            
            %%
%             qx0=zeros(MUA.Nnodes,1) ;
%             qy0=zeros(MUA.Nnodes,1) ; 
%             qx1=zeros(MUA.Nnodes,1) ;
%             qy1=zeros(MUA.Nnodes,1) ; 
            [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistent(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1);
            %%
        case "inconsistent"
            
            %%
            [UserVar,rh,kv,Tv,Lv,Pv]=LevelSetEquationAssemblyNR2inconsistent(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1);
            %%
            
            Qx=[] ; Qy=[] ; Rv=[] ; 
            
        otherwise
            
            error("case not found")
            
    end
        
    
    
   
    
end