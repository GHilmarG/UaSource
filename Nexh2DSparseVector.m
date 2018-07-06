function [h1,lambdah]=Nexh2DSparseVector(dt,h0,u0,v0,a0,u1,v1,a1,coordinates,connectivity,nip,L,Lrhs,lambdah,CtrlVar)
	
narginchk(15,15)
	
	[kv,rh]=Next2DAssembleMatrix(dt,h0,u0,v0,a0,u1,v1,a1,coordinates,connectivity,nip,CtrlVar);
	[h1,lambdah]=solveKApe(kv,L,rh,Lrhs,h0,lambdah,CtrlVar);
	h1=full(h1);

	
end




