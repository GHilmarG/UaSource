function [Tx,Fx,Ty,Fy,vs]=uvKernelFull(MUA,Deriv,DetJ,conn,sv,hv,uv,vv,AGlenv,nv,Cv,mv,qv,mukv,V0v,Bv,Sv,rhov,uov,vov,Cov,mov,uav,vav,Cav,mav,R,rhow,g,ca,sa,invPtr,invIdx,perm)
% Entry point for code generation of the uv kernel, see BuildUaCompiled.m. R holds the run-time options and parameters
% (MakeKernelParamsFlex.m). Ronly=false and the derivative flags are compile-time constants.
CV=R;
CV.EffectiveViscosity.CalculateDerivatives=true;
[Tx,Fx,Ty,Fy,vs]=uvElementKernel(MUA,Deriv,DetJ,conn,sv,hv,uv,vv,AGlenv,nv,Cv,mv,qv,mukv,V0v,Bv,Sv,rhov,uov,vov,Cov,mov,uav,vav,Cav,mav,CV,rhow,g,false,ca,sa,invPtr,invIdx,perm);
end
