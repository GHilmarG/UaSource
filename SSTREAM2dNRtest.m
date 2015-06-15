function [u,v,lambda,K,R]=...
		SSTREAM2dNRtest(s,h,u,v,coordinates,connectivity,DTxy,InteriorNodes,BoundaryNodes,nip,GF,...
		AGlen,C,L,Lb,lambda,n,m,alpha,rho,rhow,g,Itime,NLit,InfoLevel)
	
	fprintf(' initial error in satisfying Dirichlet BC %g  \n ',norm(L*[u;v]-Lb))
	tol=NLit.tol;
	InfoLevelUzawa=0;
	% Solves the SSTREAM equations in 2D (sparse and vectorized version) using Newton-Raphson iteration
	%
	
	% lambda are the Lagrange parameters used to enfore the boundary conditions
	
	% Newton-Raohson is:
	% K \Delta x_i = -R ; x_{i+1}= x_{i}+ \Delta x_i
	% R=T-F
	% T : internal nodal forces
	% F : external nodal forces
	% K : tangent matrix, where K is the directional derivative of R in the direction (Delta u, \Delta v)
	
	
	
	if any(h<0) ; error(' thickness negative ') ; end
	
	Nnodes=max(connectivity(:)); %[Nele,nod]=size(connectivity);
	
	diffVector=zeros(100,1); diffDu=1e10; fgamma=1e10 ; iteration_max=NLit.nNR; iteration_min=2 ;
	
	
	iteration=0; Progressing=1;
	while (fgamma> tol && diffDu > tol  && iteration <= iteration_max  && Progressing==1 )  || iteration < iteration_min
		
		
		iteration=iteration+1;
		
		
		[K,R,T,F]=KRTF(s,h,u,v,AGlen,n,C,m,coordinates,connectivity,nip,GF,alpha,rho,rhow,g,NLit);
		
		
		ForceResidual=norm(T(InteriorNodes)-F(InteriorNodes))/norm(F(InteriorNodes)); % Force residual at current solution
		
		
		SolMethod='Uzawa2';
		% Assuming that[u;v] already fullfill the Dirichlet BCs, I now use L [du;dv]=0 as a BCs
		
		
		[sol,dlambda]=solveKApeSymmetric(K,L,-R-L'*lambda,Lb-L*[u;v],[u;v],lambda,iteration+Itime-1,InfoLevelUzawa,SolMethod);
		
		
		
		
		du=real(sol(1:Nnodes)) ; dv=real(sol(Nnodes+1:2*Nnodes));

		
		
		% Lb must be the change needed in  [du;dv] so that [u;v]+[du;dv]=Lb
		%
		% assuming
		% L [u+gamma du;v+ gamma dv] = Lb
		% L [u;v] + gamma L [du;dv]=Lb
		% hence L [du;dv] = (Lb -L [u;v])/gamma
		% Define dLb=(Lb-L [u;v])/gamma]
		
		
		% If I start by making sure that L [u;v]=Lb, then I simply must make sure that L [du;dv]=0 after that
		
		%[sol,lambda]=solveKApeSymmetric(K,L,-R,Lb-L*[u;v],[u;v],lambda,iteration+Itime-1,InfoLevelUzawa,SolMethod);
		
		
		% calculating force residuals for this change in velocity
		
		% getting information for line search by calculating force residuals for two further step sizes
		
		
		
		% Evaluate force residual at full Newton step
		gamma=1;
		[~,T,F]=RTF(s,h,u+gamma*du,v+gamma*dv,AGlen,n,C,m,coordinates,connectivity,nip,GF,alpha,rho,rhow,g);
		fgamma=norm(T(InteriorNodes)-F(InteriorNodes))/norm(F(InteriorNodes)); % Force residuals at c times the Newton step
		
	
		
		
		if fgamma/ForceResidual >0.5
			
			
			fprintf(' Full Newton step gives a reduction of only %g so I try a line search \n',fgamma/ForceResidual)
			
			a=0 ; fa=ForceResidual;  b=0.5 ; c=gamma; fc=fgamma; solrange=[0 ;5] ; solbracket=solrange;
			
			
			[~,T1,F1]=RTF(s,h,u+b*du,v+b*dv,AGlen,n,C,m,coordinates,connectivity,nip,GF,alpha,rho,rhow,g);
			fb=norm(T1(InteriorNodes)-F1(InteriorNodes))/norm(F1(InteriorNodes)); % Force residual if full Newton Step accepted
			
			
			fvector=[fa;fb;fc]; distancevector=[a;b;c];
			
			if fa+(fc-fa)*b/c < fb  && fa < fb
				fprintf(' curvature not positive and to the right of min , backtrack! \n ')  ;
				gamma=(a+b)/2; solbracket(2)=b;
			else
				gamma = parabolamin(a,b,c,fa,fb,fc); gamma=min([max([solrange(1) gamma]) solrange(2)]);
			end
			
			
			[~,T,F]=RTF(s,h,u+gamma*du,v+gamma*dv,AGlen,n,C,m,coordinates,connectivity,nip,GF,alpha,rho,rhow,g);
			fgamma=norm(T(InteriorNodes)-F(InteriorNodes))/norm(F(InteriorNodes));
			
			fvector=[fvector ; fgamma ] ; distancevector=[distancevector ; gamma];
			[distancevector,isort]=sort(distancevector); fvector=fvector(isort) ;
			
			[smallestf,imin]=min(fvector); bestdistance=distancevector(imin);
			
			
			% gamma is the step length based on quadradic approximation
			% If the Newton-Raphson assumptions are fullfilled, we have gamma=1
			% usually this is the case exepct sometimes in the beginning and towards the end if residuals become comparable to
			% machine precision
			
			
			gammalast=-1000; tolLineSearch=0.001;
			
			ItNl=0;
			while smallestf > ForceResidual/2 || ItNl < 5  || (smallestf > ForceResidual && ItNl <10 ) % only go in hear if f is not reduced by 50% in first atempt
				ItNl=ItNl+1;
				
				if imin>1 && imin<length(fvector)
					% fprintf(' case abc \n')
					solbracket=[distancevector(imin-1) ; distancevector(imin+1)] ;
					a=distancevector(imin-1) ; fa=fvector(imin-1);
					b=distancevector(imin)   ; fb=fvector(imin);
					c=distancevector(imin+1) ; fc=fvector(imin+1);
					gamma = parabolamin(a,b,c,fa,fb,fc); gamma=min([max([solrange(1) gamma]) solrange(2)]);
					% check if new gamma to close to old value
					if abs(gamma-gammalast)< (distancevector(2)-distancevector(1))/10;
						%fprintf(' bracketing ')
						if fa > fc
							gamma=(a+b)/2;
						else
							gamma=(b+c)/2;
						end
					end
				elseif imin==1
					fprintf(' imin=1 \n')
					dd=unique(sort(distancevector));
					if dd(2)< solbracket(2)
						solbracket(2)=dd(2) ;
					end
					gamma=(dd(1)+dd(2))/2;
				elseif imin==length(fvector)
					%fprintf(' imin=length(fvector) \n')
					if distancevector(imin) < solbracket(2) && distancevector(imin) > solbracket(1) ;
						solbracket(1)=distancevector(imin) ;
					elseif distancevector(imin-1) < solbracket(2) && distancevector(imin-1) > solbracket(1) ;
						solbracket(1)=distancevector(imin-1) ;
					end
					gamma=mean(solbracket);
				end
				
				if gamma < solbracket(1) || gamma > solbracket(2) ; fprintf(' gamma put within solbracet \n') ; gamma=mean(solbracket); end
				
				gamma=min([max([solrange(1) gamma]) solrange(2)]);
				
				[~,T,F]=RTF(s,h,u+gamma*du,v+gamma*dv,AGlen,n,C,m,coordinates,connectivity,nip,GF,alpha,rho,rhow,g);
				fgamma=norm(T(InteriorNodes)-F(InteriorNodes))/norm(F(InteriorNodes));
				
				fprintf(' %g < %g < %g with fgamma % g , smallest f %g and ForceResidual %g \n ',solbracket(1),gamma,solbracket(2),fgamma,smallestf,ForceResidual)
				
				fvector=[fvector ; fgamma ] ; distancevector=[distancevector ; gamma];
				[distancevector,isort]=sort(distancevector); fvector=fvector(isort) ;
				
				
				if solbracket(2)-solbracket(1) < tolLineSearch
					
					fprintf(' tolerance in bracketing reached, breaking \n')
					
					
					break
				end
				if abs(gamma-gammalast)< tolLineSearch && smallestf < ForceResidual 
					fprintf(' change in gamma too small , breaking \n')
					break
				end
				
				if ItNl > 2 && gamma~=0 && smallestf < ForceResidual 
					smallestvalues=unique(sort(fvector));
					if (smallestvalues(2)-smallestvalues(1))/ForceResidual <0.01
						fprintf(' further reduction in function too small, breaking after %g iterations \n',ItNl)
						break
					end
				end
				
				[smallestf,imin]=min(fvector); bestdistance=distancevector(imin);
				gammalast=gamma;
				
				
			end
			
			
			
			if fgamma>smallestf && bestdistance>=solbracket(1) && bestdistance <= solbracket(2)
				gamma=bestdistance ; fgamma=smallestf ;
			end
			
			if gamma==0 ;
				fprintf(' gamma returned is zero ! Algorithim has stagnated \n ') ;
				Progressing=0;
			end
			
			[~,T,F]=RTF(s,h,u+gamma*du,v+gamma*dv,AGlen,n,C,m,coordinates,connectivity,nip,GF,alpha,rho,rhow,g);
			fgamma=norm(T(InteriorNodes)-F(InteriorNodes))/norm(F(InteriorNodes)); % Force residuals at gamma times the Newton step
			
			figure(999) ; hold off ; plot(distancevector,fvector,'o') ; hold on ; plot(gamma,fgamma,'+r') ; plot(solbracket,[0;0],'xr') ;
			
			
			%        % step size selected by requiring residual force at the end of iteration to be orthogonal to (\Delta u, \Delta v)
			% 		R0=[du ; dv]'*R; R1=[du ; dv]'*R1;
			% 		a=R0/R1; asqr=sqrt((a/2)^2-a); ap=a/2+asqr;  am=a/2-asqr;
			% 		if a<0 ; gamma=ap ; else gamma=a/2; end
			% 		gamma=max([0.25 min([gamma 1])]);  % just in case quadradic appoximation is totally off
			% 		fprintf(' gamma %g R0 %g R1 %g ap %g am %g \n',gamma,R0,R1,ap,am)
			
			if InfoLevel > 0 
				fprintf(' Line search gives gamma % g and Force residual for [u;v]+gamma [du;dv] %g , compared to %g for zero step size \n ',...
					gamma,fgamma, ForceResidual)
				
			end
		end
		
		
		
		D=mean(sqrt(u.*u+v.*v));
		diffDu=(max(abs(gamma*du))+max(abs(gamma*dv)))/D; % sum of max change in du and dv normalized by mean speed
		
		
		diffVector(iteration)=fgamma;
		%diffVector(iteration)=diffDu;
		
		
		
		u=u+gamma*du ; v=v+gamma*dv;
		lambda=lambda+gamma*dlambda;
		
		figure
	    PlotForceResidualVectors(T-F,coordinates,InteriorNodes)
		
		
		fprintf('NR iteration # %g , Last ForceResidual %g , New ForceResidual %g , Ratio %g, diffDu %g , gamma %g \n ',...
			iteration, ForceResidual,fgamma,fgamma/ForceResidual,diffDu,gamma)
		
		
		
	end
	
	if InfoLevel>5 && iteration >= 2
		
		
		N=max([1,iteration-5]);
		
		[detrended,a0,a1]=detrend_xt(log10(diffVector(N:iteration)),N:iteration);
		fprintf(' slope NR : %g \n',a1)
		figure(2100) ;
		semilogy(diffVector(1:iteration),'x-g') ; title('NR') ;
		% hold on ; semilogy(10.^(a0+a1*[1:iteration]),'ko') ;
		
	end
	
	if iteration > iteration_max
		warning('SSTREAM2dNR:MaxIterationReached','SSTREAM2NR exits because maximum number of iterations %g reached \n',iteration_max)
	end
	
end

