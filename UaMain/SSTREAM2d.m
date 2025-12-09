

function [u,v,lambda,kv,rh]=...
		SSTREAM2d(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,C,L,Lb,lambda,n,m,alpha,rho,rhow,g,CtrlVar,Itime)
	
	
	

		[u,v,lambda,kv,rh]=...
			SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,C,L,Lb,lambda,n,m,alpha,rho,rhow,g,Itime,CtrlVar);

	
	
end

