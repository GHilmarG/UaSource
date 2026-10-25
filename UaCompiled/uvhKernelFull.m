function [Tx,Fx,Ty,Fy,Th,Fh,Xmat,tauxInt,tauyInt,etaInt,HeInt]=uvhKernelFull(MUA,conn,bv,hv,uv,vv,AGlenv,nv,Cv,mv,qv,mukv,V0v,h0v,u0v,v0v,as0v,ab0v,as1v,ab1v,dadhv,Bv,Sv,rhov,LSFMaskv,hBCsMaskv,uov,vov,Cov,mov,uav,vav,Cav,mav,R,rhow,g,ca,sa,dt)
% Entry point for code generation, see BuildUaCompiled.m. R holds the run-time options and parameters (MakeKernelParamsFlex.m).
CV=R;
CV.BasalDrag.CalculateDerivatives=true;
CV.EffectiveViscosity.CalculateDerivatives=true;
[Tx,Fx,Ty,Fy,Th,Fh,Xmat,tauxInt,tauyInt,etaInt,HeInt]=...
    uvhElementKernel(MUA,conn,bv,hv,uv,vv,AGlenv,nv,Cv,mv,qv,mukv,V0v,h0v,u0v,v0v,as0v,ab0v,as1v,ab1v,dadhv,Bv,Sv,rhov,LSFMaskv,hBCsMaskv,uov,vov,Cov,mov,uav,vav,Cav,mav,CV,rhow,g,false,ca,sa,dt);
end
