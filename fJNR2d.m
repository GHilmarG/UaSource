function [J,u,v]=fJNR2d(gamma,dJdC,dJdA,uMeas,vMeas,wMeasInt,s,S,B,h,u,v,wint,Cd,CAGlen,CC,AGlen_prior,C_prior,coordinates,connectivity,Boundary,nip,...
		AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar)
	
	persistent ulast vlast
    
	% Returns value of Misfit function J as a function of gamma by first calculating u and v
    % positivity in C and A ensured by using the projected gradient method
	
	
	if ~isempty(ulast)
		u=ulast ;v=vlast;
		%fprintf(' fJNR2d using locally persistent values for u and v.   ')
	end
	
	[C]=ProjGradient(C,dJdC,gamma,CtrlVar.Cmin,CtrlVar.Cmax);
	[AGlen]=ProjGradient(AGlen,dJdA,gamma,CtrlVar.AGlenmin,CtrlVar.AGlenmax);
	
	fprintf(' fJNR2d solving forward problem \n ')
	[u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
	
	ulast=u; vlast=v ;
	
    J=MisfitFunction(u,v,wint,uMeas,vMeas,wMeasInt,C,C_prior,AGlen,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,CtrlVar)    ;
    
    
	
end

