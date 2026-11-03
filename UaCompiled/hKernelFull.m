function [f0,vs]=hKernelFull(MUA,Deriv,DetJ,conn,h0v,h1v,u0v,v0v,u1v,v1v,rhov,a0v,a1v,dadhv,LSFMaskv,hBCsMaskv,R,rPtr,rIdx,invPtr,invIdx)
% Entry point for code generation of the h kernel, see BuildUaCompiled.m. R holds the run-time options and parameters
% (MakeKernelParamsFlex.m). Ronly=false is a compile-time constant.
[f0,vs]=hElementKernel(MUA,Deriv,DetJ,conn,h0v,h1v,u0v,v0v,u1v,v1v,rhov,a0v,a1v,dadhv,LSFMaskv,hBCsMaskv,R,false,rPtr,rIdx,invPtr,invIdx);
end
