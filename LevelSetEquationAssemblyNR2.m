function [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1)
    
  
    narginchk(15,15)
    nargoutchk(2,9) 
    
    if ~(nargout==2 || nargout==9)
        error('af')
    end
    
    
    
    switch lower(CtrlVar.LevelSetAssembly)
        
        case "consistent"
            
            
            if CtrlVar.Parallel.LSFAssembly.parfor.isOn
                [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistentParfor(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1);
            else
                [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistent(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1);
            end
            
             CtrlVar.Parallel.isTest=true; 
            if CtrlVar.Parallel.isTest
                tic;
                [UserVar,rhT,kvT,TvT,LvT,PvT,QxT,QyT,RvT]=LevelSetEquationAssemblyNR2consistentParfor(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1);
                tSeq=toc;
                
                tic;
                [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy,Rv]=LevelSetEquationAssemblyNR2consistent(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1,qx0,qy0,qx1,qy1);
                tPar=toc;
               
                fprintf("Testing: norm(full(rh-rhT))=%g \n ",norm(full(rh-rhT))/norm(rh))
                fprintf("Testing: norm(full(diag(kv)-diag(kvT)))=%g \n ",norm(diag(kv)-diag(kvT))/norm(diag(kv)))
                
                fprintf("\n parfor \t\t\t  for \t\t\t speedup \n")
                fprintf("%f \t\t\t %f \t\t %f \n",tPar,tSeq,tSeq/tPar)
                
            end
            
        case "inconsistent"
            
            %%
            [UserVar,rh,kv,Tv,Lv,Pv]=LevelSetEquationAssemblyNR2inconsistent(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1);
            %%
            
            Qx=[] ; Qy=[] ; Rv=[] ; 
            
        otherwise
            
            error("case not found")
            
    end
        
    
    
   
    
end