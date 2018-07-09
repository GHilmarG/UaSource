function [beta2,Dbeta2Duu,Dbeta2Dvv,Dbeta2Duv] = calcBeta2in2D(u,v,C,m,GF)
	
	global SpeedZero beta2Max beta2Min
	
		
	% calculated beta^2 and D beta^2/Du at nodal points
	
	

	beta2= C.^(-1/m).*(sqrt(u.*u+v.*v+SpeedZero^2)).^(1/m-1) ;
	
	beta2(beta2>beta2Max)=beta2Max; beta2(beta2<beta2Min)=beta2Min;
	
	beta2=GF.node.*beta2;
	
	if nargout>1
		
		% The directional derivative is
		% D beta^2(u,v)[Delta u, Delta v]= (1/m-1) C^(-1/m) (u^2+v^2)^((1-3m)/2m)  (u \Delta u + v \Delta v)
		%
		
		Dbeta2=(1/m-1).*C.^(-1/m).*(u.^2+v.^2+SpeedZero^2).^((1-3*m)/(2*m));
		
		Dbeta2=Dbeta2.*GF.node;
		
		Dbeta2Duu=Dbeta2.*u.*u;
		Dbeta2Dvv=Dbeta2.*v.*v;
		Dbeta2Duv=Dbeta2.*u.*v;

	end
	
	
end

