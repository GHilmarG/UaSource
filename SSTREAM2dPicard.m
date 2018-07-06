

function [u,v,lambda,kv,rh,etaInt]=...
		SSTREAM2dPicard(s,S,h,u,v,coordinates,connectivity,Boundary,nip,etaInt,GF,AGlen,C,L,Lb,lambda,n,m,alpha,rho,rhow,g,Itime,NLit,InfoLevel)
                      
    
    % Solves the SSTREAM equations in 2D (sparse and vectorized version)
    
    % lambda are the Lagrange parameters used to enfore the boundary conditions
    
    if any(h<0) ; error(' thickness negative ') ; end
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    
    diff=1e10;
    
        
    
    if n==1 && m==1 ;
        iteration_max=1 ; iteration_min=1 ; nonlinear=0;
    else
        fprintf(' non-linear n %g m %g \n',n,m)
        iteration_max=NLit.nP; iteration_min=2 ; nonlinear=1;
    end
    
    iteration=0;
    while (diff > NLit.tol  && iteration < iteration_max )  || iteration < iteration_min
        iteration=iteration+1;
        %disp(' ') ; disp('-------------------------------------------------------------------------')
        %disp([' SSS iteration # : ',num2str(iteration)])
        
		beta2 = calcBeta2in2D(u,v,C,m,GF);
        %beta2= C.^(-1/m).* (sqrt(u.*u+v.*v)).^(1/m-1) ;
        		
		tassemble=tic;
        [kv,rh]=kvrhPicard(s,h,coordinates,connectivity,nip,etaInt,beta2,alpha,rho,rhow,g);
        tassemble=toc(tassemble);
		
        tStartSolve=tic;
		SolMethod='Uzawa2';
        [sol,lambda]=solveKApeSymmetric(kv,L,rh,Lb,[u;v],lambda,iteration+Itime-1,InfoLevel,SolMethod);
        tElapsedSolve=toc(tStartSolve);
		
      
        
        ulast=u ; vlast=v ;  
		%u=NLit.Pgamma*real(sol(1:Nnodes))+(1-NLit.Pgamma)*ulast ; v=NLit.Pgamma*real(sol(Nnodes+1:2*Nnodes))+(1-NLit.Pgamma)*vlast;
		u=real(sol(1:Nnodes)) ; v=real(sol(Nnodes+1:2*Nnodes));
		
        
        if nonlinear
            D=mean(sqrt(u.*u+v.*v));
            diffv=norm(v-vlast)/sqrt(length(v));  % rms of change in u and v
            diffu=norm(u-ulast)/sqrt(length(u));
            diff=diffu/D+diffv/D;                 % normalized by the mean speed
            
            if InfoLevel > 0 ;fprintf('Picard iteration # %g , diff %g , soltime %g , assemble time %g : \n ',iteration, diff,tElapsedSolve,tassemble) ; end
            
            if n~=1
				etaIntLast=etaInt;
                [etaInt]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n);
				etaInt=NLit.Pgamma*etaInt+(1-NLit.Pgamma)*etaIntLast;
            end
            
        end
    end
end

