function [h1,lambdah]=SSHEETtransient2HD(CtrlVar,connectivity,coordinates,AGlen,n,C,m,rho,rhow,g,s0,b0,B,S,u0,v0,a0,a1,nip,L,Lrhs,lambdah)
	

	
	[kv,rh]=MatrixAssemblySSHEETtransient2HD(CtrlVar,connectivity,coordinates,AGlen,n,C,m,rho,rhow,g,s0,b0,B,S,u0,v0,a0,a1,nip,L,Lrhs,lambdah);
	[h1,lambdah]=solveKApe(kv,L,rh,Lrhs,h0,lambdah,CtrlVar);
	h1=full(h1);

	
end
