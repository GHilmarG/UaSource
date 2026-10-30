function [Tint,Fext,vs,tauxInt,tauyInt,etaInt,HeInt]=uvhKernelRonly(MUA,Deriv,DetJ,conn,bv,hv,uv,vv,AGlenv,nv,Cv,mv,qv,mukv,V0v,h0v,u0v,v0v,as0v,ab0v,as1v,ab1v,dadhv,Bv,Sv,rhov,LSFMaskv,hBCsMaskv,uov,vov,Cov,mov,uav,vav,Cav,mav,R,rhow,g,ca,sa,dt,IntOut,rPtr,rIdx,invPtr,invIdx)
% Entry point for code generation, see BuildUaCompiled.m. R holds the run-time options and parameters (MakeKernelParamsFlex.m).
% Ronly=true and the derivative flags are compile-time constants. (invPtr and invIdx are not used.)
CV=R;
CV.EffectiveViscosity.CalculateDerivatives=false;
[Tint,Fext,vs,tauxInt,tauyInt,etaInt,HeInt]=...
    uvhElementKernel(MUA,Deriv,DetJ,conn,bv,hv,uv,vv,AGlenv,nv,Cv,mv,qv,mukv,V0v,h0v,u0v,v0v,as0v,ab0v,as1v,ab1v,dadhv,Bv,Sv,rhov,LSFMaskv,hBCsMaskv,uov,vov,Cov,mov,uav,vav,Cav,mav,CV,rhow,g,true,ca,sa,dt,IntOut,rPtr,rIdx,invPtr,invIdx);
end
