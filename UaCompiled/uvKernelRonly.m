function [Tx,Fx,Ty,Fy,vs]=uvKernelRonly(MUA,Deriv,DetJ,conn,sv,hv,uv,vv,AGlenv,nv,Cv,mv,qv,mukv,V0v,Bv,Sv,rhov,uov,vov,Cov,mov,uav,vav,Cav,mav,R,rhow,g,ca,sa,invPtr,invIdx,perm)
% Entry point for code generation of the uv kernel, see BuildUaCompiled.m. R holds the run-time options and parameters
% (MakeKernelParamsFlex.m). Ronly=true and the derivative flags are compile-time constants. (invPtr, invIdx and perm are not used.)
CV=R;
CV.EffectiveViscosity.CalculateDerivatives=false;
[Tx,Fx,Ty,Fy,vs]=uvElementKernel(MUA,Deriv,DetJ,conn,sv,hv,uv,vv,AGlenv,nv,Cv,mv,qv,mukv,V0v,Bv,Sv,rhov,uov,vov,Cov,mov,uav,vav,Cav,mav,CV,rhow,g,true,ca,sa,invPtr,invIdx,perm);
end
