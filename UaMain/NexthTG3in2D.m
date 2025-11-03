function [h1,lambdah]=NexthTG3in2D(dt,h0,u0,v0,du0dt,dv0dt,a0,da0dt,u1,v1,a1,da1dt,du1dt,dv1dt,coordinates,connectivity,Boundary,nip,L,Lrhs,lambdah,CtrlVar)
	
% Third order impicit Taylor-Galerkin (implict with respect to h, not u and v)
	

	[kv,rh]=NexthTG3Assemble2DMatrix(dt,h0,u0,v0,du0dt,dv0dt,a0,da0dt,u1,v1,a1,da1dt,du1dt,dv1dt,coordinates,connectivity,Boundary,nip,CtrlVar);
	
	[h1,lambdah]=solveKApe(kv,L,rh,Lrhs,h0,lambdah,CtrlVar);
	h1=full(h1);

	
end
