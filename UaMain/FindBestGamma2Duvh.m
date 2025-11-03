function [r,gamma,ExitFlag] = FindBestGamma2Duvh(F0,raccept,r0,r1,u,v,h,du,dv,dh,S,B,u0,v0,h0,L,lambda,dlambda,as,ab,dt,AGlen,n,C,m,coordinates,connectivity,nip,alpha,rho,rhow,g,MeshProp,CtrlVar)
	
	
    %
    % ExitFlag : 0 if everything OK
    %            1 if final value not better than starting value
    
    ExitFlag=0;
    
	Slope=-2*r0 ; r=r1 ;
	
	
	%   Assume the full Newton step has already been tried and rejected
	%   Always try backtracking, and if that does not give acceptable fit try general line search
	%
	%if (r >  r0-beta*Slope*gamma)  % if not sufficient decrease try cubic backtracking
	
	
    
	if CtrlVar.InfoLevelNonLinIt>4
		if r >  raccept
			fprintf(' Newton step overshot with r1/r0=%-14.7g so try backtracking \n',r1/r0)
		else
			fprintf(' Full Newton step was rejected for r1/r0=%-14.7g so try backtracking \n',r1/r0)
		end
	end
	
	gamma=-Slope/2/(r1-r0-Slope);  % location of minimum based on a quadradic fit
	gamma=max([0.1 min([0.8 gamma])]);
	
	% calculate cost function at minimum for quadradic fit
	r=CalcCostFunctionNRuvh(gamma,du,dv,dh,u,v,h,S,B,u0,v0,h0,as,ab,dt,AGlen,n,C,m,coordinates,connectivity,nip,alpha,rho,rhow,g,MeshProp,CtrlVar,F0,L,lambda,dlambda);
	
	
	if CtrlVar.InfoLevelNonLinIt>4
		fprintf('                                 After quadradic backtracking : gamma %-14.7g r %-14.7g ratio %-14.7g \n ',gamma,r,r/r0)
		
	end
	
	% cubic backtracking only requires one additional evaluation of r, so always do that as well and the
	% select the best value from quadradic and cubic backtracking
	%if r > ( r0-beta*Slope*gamma) || (r/r0> raccept)
	% Try cubic fit
	if CtrlVar.InfoLevelNonLinIt>4 ; fprintf(' Try cubic backtracking \n' ) ; end
	rb=r ; b=gamma;
	% and find minimum based on a cubic fit
	[gamma] = CubicFit(Slope,r0,rb,r1,b,1); gamma=max([0.1*b min([0.5*b gamma])]);
	r=CalcCostFunctionNRuvh(gamma,du,dv,dh,u,v,h,S,B,u0,v0,h0,as,ab,dt,AGlen,n,C,m,coordinates,connectivity,nip,alpha,rho,rhow,g,MeshProp,CtrlVar,F0,L,lambda,dlambda);
	
	
	if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('Cubic backtracking fit gives : gamma %-14.7g r %-14.7g ratio %-14.7g \n',gamma,r,r/r0) ; end
	
	if r > rb
		if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('Cubic backtracking fit gives worse results than quadractic fit, i.e. %-g instead of %-g with gamma %g\n',r,rb,gamma) ; end
		gamma=b ; r=rb ; % just use values from quadradic fit
	else
		if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('Cubic backtracking fit gives better results than quadractic fit, i.e. %-g instead of %-g with gamma %g\n',r,rb,gamma) ; end
	end
	
	if  r >  raccept  % if the ratio is still too small, I might have a problem, rather than giving up, try general line search
		
		if CtrlVar.InfoLevelNonLinIt>3
			fprintf(' try line search as r/r0=%-g is still too big despite having tried backtracking \n ',r/r0)
		end
		
		
		[gamma,r]=fminbndGHG(@(gamma) ...
			CalcCostFunctionNRuvh(gamma,du,dv,dh,u,v,h,S,B,u0,v0,h0,as,ab,dt,AGlen,n,C,m,coordinates,connectivity,nip,alpha,rho,rhow,g,MeshProp,CtrlVar,F0,L,lambda,dlambda),...
			CtrlVar.NewtonMinStep,CtrlVar.NewtonMaxStep,1,r1,optimset('TolX',0.001,'Display','off','MaxFunEvals',10));
		
		
		if CtrlVar.InfoLevelNonLinIt>4
			fprintf('After line search: gamma=%-14.7g r=%-14.7g, ratio=%-14.7g \n ',gamma,r,r/r0)
			
		end
	end
	
	
	if r>r0 ;
		fprintf('FindBestGamma:  best value not better than starting value %-g to %-g \n !',r0,r)
        ExitFlag=1;
    end
	
end

