function [f0,vs]=hKernelRonly(MUA,Deriv,DetJ,conn,h0v,h1v,u0v,v0v,u1v,v1v,rhov,a0v,a1v,dadhv,LSFMaskv,hBCsMaskv,R,rPtr,rIdx,invPtr,invIdx)
% Entry point for code generation of the h kernel, see BuildUaCompiled.m. R holds the run-time options and parameters
% (MakeKernelParamsFlex.m). Ronly=true is a compile-time constant. (invPtr and invIdx are not used.)
[f0,vs]=hElementKernel(MUA,Deriv,DetJ,conn,h0v,h1v,u0v,v0v,u1v,v1v,rhov,a0v,a1v,dadhv,LSFMaskv,hBCsMaskv,R,true,rPtr,rIdx,invPtr,invIdx);
end
