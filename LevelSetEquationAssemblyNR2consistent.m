function [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistent(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1)
    
 
    narginchk(15,15)
    
    CtrlVar.LevelSetMethodSourceTermFomulation=true; 
    if CtrlVar.LevelSetMethodSourceTermFomulation
        
          [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistentST(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1);
        
    else
        
        [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistentNST(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1);
        
    end
    
    
   
end