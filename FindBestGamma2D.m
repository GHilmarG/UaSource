function [r,gamma] = FindBestGamma2D(F0,r0,r1,s,S,B,h,u,du,v,dv,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,L,lambda,dlambda,CtrlVar)
	
	beta=1e-4;  % Only accept step if reduction is larger than beta*StepSize
	Slope=-2*r0 ;  % using the inner product def
	
	
	r=r1 ; gamma=1;
	
	if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('\n \n \n') ; end
	
    target= r0-beta*Slope*gamma  ;
    
    if r >  target || (r/r0>CtrlVar.NewtonAcceptRatio)  % if not sufficient decrease try  backtracking
        
        if CtrlVar.InfoLevelNonLinIt>4
            fprintf(' Full Newton step with r1/r0=%-14.7g not accepted, start backtracking \n',r1/r0)
        end
        
        % initial quadric backtracking
        gammab=-Slope/2/(r1-r0-Slope); gammab=max([0.1 min([0.8 gammab])]);
        rb = CalcCostFunctionNR(gammab,s,S,B,h,u,du,v,dv,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,F0,L,lambda,dlambda,CtrlVar);
        
        if CtrlVar.InfoLevelNonLinIt>4
            fprintf('                                 After quadradic backtracking : \t gamma=%-14.7g \t r=%-14.7g \t ratio=%-14.7g \n ',gammab,rb,rb/r0)
        end
        
        if rb > r
            if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('quadradic backtracking fit rejected \n'); end
            
            % if quadracic backtracking rejected, try parabolic fit through these three functional values:
            gammap=parabolamin(0,gammab,1,r0,rb,r1);
            rp = CalcCostFunctionNR(gammap,s,S,B,h,u,du,v,dv,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,F0,L,lambda,dlambda,CtrlVar);
            if CtrlVar.InfoLevelNonLinIt>4 ;
                fprintf('                                 After quadradic fit : \t gamma=%-14.7g \t r=%-14.7g \t ratio=%-14.7g \n ',gammap,rp,rp/r0)
            end
            if rp<r
                r=rp ; gamma=gammap;
                if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('quadradic  fit accepted \n'); end
                
                if gammap>1 % if accepted and to the right of Newton step, try one additional quadradtic fit.
                    gammap2=parabolamin(gammab,1,gammap,rb,r1,rp);
                    rp2 = CalcCostFunctionNR(gammap2,s,S,B,h,u,du,v,dv,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,F0,L,lambda,dlambda,CtrlVar);
                    if CtrlVar.InfoLevelNonLinIt>4 ;
                        fprintf('                       After a further quadradic fit : \t gamma=%-14.7g \t r=%-14.7g \t ratio=%-14.7g \n ',gammap2,rp2,rp2/r0)
                    end
                    if rp2<r
                        r=rp2 ; gamma=gammap2;
                        
                        if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('a further quadradic  fit accepted \n'); end
                    end
                end
            end
        else
            r=rb ; gamma=gammab;
            if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('quadradic backtracking fit accepted \n'); end
            target= r0-beta*Slope*gamma  ;
            
            if r >  target ||( r/r0>CtrlVar.NewtonAcceptRatio)
                % if quadradic fit accepted, but still not sufficiten decreas, try cubic backtracking
                if CtrlVar.InfoLevelNonLinIt>4 ; fprintf(' Try cubic backtracking \n' ) ; end
                
                % and find minimum based on a cubic fit
                [gammac] = CubicFit(Slope,r0,rb,r1,gammab,1); gammac=max([0.1*gammab min([0.5*gammab gammac])]);
                rc = CalcCostFunctionNR(gammac,s,S,B,h,u,du,v,dv,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,F0,L,lambda,dlambda,CtrlVar);
                
                if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('After cubic backtracking : \t gamma=%-14.7g \t r=%-14.7g \t ratio=%-14.7g \n',gammac,rc,rc/r0) ; end
                
                if rc > r
                    if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('Cubic backtracking fit rejected \n'); end
                else
                    r=rc ; gamma=gammac;
                    if CtrlVar.InfoLevelNonLinIt>4 ; fprintf('Cubic backtracking fit accepted \n'); end
                end
            end
        end
        
        
    end
	
    target= r0-beta*Slope*gamma  ;
	if r > target % still not sufficient decrease despite having tried cubic backtracking, line search
		
		if CtrlVar.InfoLevelNonLinIt>3
			fprintf(' try line search as ratio is %g \n ',r/r0)
		end
		gammaOld=gamma ; rOld=r;
		
		[gamma,r]=fminbndGHG(@(gamma) CalcCostFunctionNR(gamma,s,S,B,h,u,du,v,dv,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,F0,L,lambda,dlambda,CtrlVar),...
			CtrlVar.NewtonMinStep,CtrlVar.NewtonMaxStep,gamma,r,optimset('TolX',CtrlVar.NewtonMinStep,'Display','off','MaxFunEvals',20));
		
		
		if CtrlVar.InfoLevelNonLinIt>4
			fprintf('After line search: \t \t \t gamma=%-14.7g \t r=%-14.7g \t ratio=%-14.7g \n ',gamma,r,r/r0)
			
		end
		
		if r> rOld ; r=rOld ; gamma=gammaOld ; end
	end
	
	
	if r>r0 && r>CtrlVar.NLtol; % if everyting else fails, allow some increase, possibly it must get out of a local minimum
		fprintf('FindBestGamma: residual increased from %g to %g, but still returning the new value \n !',r0,r)
		% 		if r>10*r0
		% 			fprintf('FindBestGamma: residual increased from %g to %g, returning original values \n !',r0,r)
		% 			r=r0; gamma=0;
		% 		end
	end
	
end

