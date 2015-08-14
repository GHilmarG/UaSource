function  [ub,vb,ubvbLambda,Kuv,Ruv,RunInfo,ubvbL]=SSTREAM2dNR(CtrlVar,MUA,BCs,s,S,B,h,ub,vb,ubvbLambda,AGlen,C,n,m,alpha,rho,rhow,g)
                   
    narginchk(18,18)
 
    tStart=tic;
    RunInfo.converged=1; RunInfo.Iterations=NaN;  RunInfo.residual=NaN;
    
    MLC=BCs2MLC(MUA,BCs);
    ubvbL=MLC.ubvbL; ubvbRhs=MLC.ubvbRhs;

    
    
    
    if isempty(ubvbRhs)
        ubvbLambda=[];
    elseif numel(ubvbLambda)~=numel(ubvbRhs) ;
        ubvbLambda=zeros(numel(ubvbRhs),1) ;
    end
    
    if any(isnan(C)) ; save TestSave ; error( ' C nan ') ; end
    if any(isnan(AGlen)) ; save TestSave ;error( ' AGlen nan ') ; end
    if any(isnan(S)) ; save TestSave ; error( ' S nan ') ; end
    if any(isnan(h)) ; save TestSave error( ' h nan ') ; end
    if any(isnan(ub)) ; save TestSave ; error( ' ub nan ') ; end
    if any(isnan(vb)) ; save TestSave ; error( ' vb nan ') ; end
    if any(isnan(ubvbLambda)) ; save TestSave ; error( ' ubvbLambda nan ') ; end
    if any(isnan(rho)) ; save TestSave  ; error( ' rho nan ') ; end
     if any(h<0) ; warning('MATLAB:SSTREAM2dNR:hnegative',' thickness negative ') ; end
    
     
 
    %%
    
    
    % Solves the SSTREAM equations in 2D (sparse and vectorized version) using Newton-Raphson iteration
    
    % lambda are the Lagrange parameters used to enfore the boundary conditions
    
    % Newton-Raphson is:
    % K \Delta x_i = -R ; x_{i+1}= x_{i}+ \Delta x_i
    % R=T-F
    % T : internal nodal forces
    % F : external nodal forces
    % K : tangent matrix, where K is the directional derivative of R in the direction (Delta u, \Delta v)
    
    % I need to solve
    %
    % [Kxu Kxv Luv'] [du]        =  [ -Ru ] - Luv' lambdauv
    % [Kyu Kyv     ] [dv]        =  [ -Rv ]
    % [  Luv      0] [dlambdauv]    [ Lrhsuv-Luv [u ;v ]
    %
    % All matrices are Nnodes x Nnodes, apart from:
    % Luv is #uv constraints x 2 Nnodes
    %
    
    
    
    
    
    
    diffVector=zeros(CtrlVar.NRitmax+1,1); diffDu=1e10; r=1e10 ;
    
    dub=zeros(MUA.Nnodes,1) ; dvb=zeros(MUA.Nnodes,1) ; dlambda=zeros(numel(ubvbLambda),1);
    
    
    
    iteration=0;
    while ((r> CtrlVar.NLtol  || diffDu > CtrlVar.du  )&& iteration <= CtrlVar.NRitmax  )  || (iteration < CtrlVar.NRitmin && (n~=1 || m~=1))
  
        
        iteration=iteration+1;
        
        if CtrlVar.CalvingFrontFullyFloating
            [Kuv,Ruv,~,F]=KRTF(s,S,B,h,ub,vb,AGlen,n,C,m,MUA.coordinates,MUA.connectivity,MUA.Boundary,MUA.nip,alpha,rho,rhow,g,CtrlVar);
        else
            [Ruv,Kuv,~,F]=KRTFgeneralBCs(CtrlVar,MUA,s,S,B,h,ub,vb,AGlen,n,C,m,alpha,rho,rhow,g);
        end
        
        F0=F;
        
        %
        % 		[Ktest,Rtest]=KRTFloop(s,S,B,h,u,v,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,CtrlVar);
        % 		save TestSave Ktest Rtest K R
        % 		error('sfda')
        
        
        if numel(ubvbL)==0
            r0=ResidualCostFunction(Ruv,F0);
        else
            r0=ResidualCostFunction(Ruv+ubvbL'*ubvbLambda,F0);
        end
        
        if ~isreal(Kuv) ; save TestSave Kuv ; error('SSTREAM2dNR: K not real') ;  end
        if ~isreal(ubvbL) ; save TestSave ubvbL ; error('SSTREAM2dNR: L not real') ;  end
        if any(isnan(Kuv)) ; save TestSave Kuv ; error('SSTREAM2dNR: K nan') ;  end
        if any(isnan(ubvbL)) ; save TestSave ubvbL ; error('SSTREAM2dNR: L nan') ;  end
        
        
        if CtrlVar.IncludeDirichletBoundaryIntegralDiagnostic==0;
            
            CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=1;
            if numel(ubvbL)==0
                [sol,dlambda]=solveKApeSymmetric(Kuv,ubvbL,-Ruv,ubvbRhs,[dub;dvb],dlambda,CtrlVar);
            else
                [sol,dlambda]=solveKApeSymmetric(Kuv,ubvbL,-Ruv-ubvbL'*ubvbLambda,ubvbRhs-ubvbL*[ub;vb],[dub;dvb],dlambda,CtrlVar);
            end

        else
            CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=0;
            [sol,dlambda]=solveKApe(Kuv,ubvbL,-Ruv-ubvbL'*ubvbLambda,ubvbRhs-ubvbL*[ub;vb],[dub;dvb],dlambda,CtrlVar);
        end
        
        

        dub=real(sol(1:MUA.Nnodes)) ; dvb=real(sol(MUA.Nnodes+1:2*MUA.Nnodes));
        
        % Evaluate force residual at full Newton step
        r1 = CalcCostFunctionNR(1,s,S,B,h,ub,dub,vb,dvb,AGlen,n,C,m,MUA,alpha,rho,rhow,g,F0,ubvbL,ubvbLambda,dlambda,CtrlVar);
        [r,gamma,infovector,BacktrackingInfo] = FindBestGamma2Dbacktracking(F0,r0,r1,s,S,B,h,ub,dub,vb,dvb,AGlen,n,C,m,MUA,alpha,rho,rhow,g,ubvbL,ubvbLambda,dlambda,CtrlVar);
        
        if BacktrackingInfo.Converged==0;    
            fprintf(CtrlVar.fidlog,' SSTREAM2dNR backtracking step did not converge \n ') ;
            warning('SSTREAM2NR:didnotconverge',' SSTREAM2dNR backtracking step did not converge \n ')
            fprintf(CtrlVar.fidlog,' saving variables in SSTREAM2dNRDump \n ') ;
            save SSTREAM2dNRDump
           
        end
        
        %% If requested, plot residual as function of steplength
        if CtrlVar.InfoLevelNonLinIt>=10 && CtrlVar.doplots==1
            nnn=10;
            gammaTestVector=zeros(nnn,1) ; rTestvector=zeros(nnn,1);
            
            Up=2.2;
            if gamma>0.7*Up ; Up=2*gamma; end
            parfor I=1:nnn
                gammaTest=Up*(I-1)/(nnn-1)+gamma/25;
                rTest=CalcCostFunctionNR(gammaTest,s,S,B,h,ub,dub,vb,dvb,AGlen,n,C,m,MUA,alpha,rho,rhow,g,F0,ubvbL,ubvbLambda,dlambda,CtrlVar);
                gammaTestVector(I)=gammaTest ; rTestvector(I)=rTest;
            end
            
            gammaTestVector=[gammaTestVector(:);infovector(:,1)];
            rTestvector=[rTestvector(:);infovector(:,2)];
            [gammaTestVector,ind]=unique(gammaTestVector) ; rTestvector=rTestvector(ind) ;
            [gammaTestVector,ind]=sort(gammaTestVector) ; rTestvector=rTestvector(ind) ;
            
            
            figure ; plot(gammaTestVector,rTestvector,'o-r') ; hold on ; 
            
            plot(gamma,r,'Marker','h','MarkerEdgeColor','k','MarkerFaceColor','g')
            
            slope=-2*r0;
            plot([gammaTestVector(1) gammaTestVector(2)],[rTestvector(1) rTestvector(1)+(gammaTestVector(2)-gammaTestVector(1))*slope],'g')
            title(sprintf('iteration %i ',iteration)) ; hold off
            
            %input('press return to continue')
        end
        
        D=mean(sqrt(ub.*ub+vb.*vb));
        diffDu=full(max(abs(gamma*dub))+max(abs(gamma*dvb)))/D; % sum of max change in du and dv normalized by mean speed
        diffDlambda=full(max(abs(gamma*dlambda))/mean(abs(ubvbLambda)));
        diffVector(iteration)=r0;   % override last value, because it was just a (very accurate) estimate
        diffVector(iteration+1)=r;
        
        ub=ub+gamma*dub ; 
        vb=vb+gamma*dvb; 
        ubvbLambda=ubvbLambda+gamma*dlambda;
        
        
        
        if CtrlVar.InfoLevelNonLinIt>100  && CtrlVar.doplots==1
            PlotForceResidualVectors('uv',Ruv,ubvbL,ubvbLambda,MUA.coordinates,CtrlVar) ; axis equal tight
        end
        if CtrlVar.InfoLevelNonLinIt>=1
            fprintf(CtrlVar.fidlog,'NRuv:%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , du=%-14.7g , dl=%-14.7g \n ',...
                iteration,BacktrackingInfo.iarm,gamma,r/r0,r0,r,diffDu,diffDlambda);
        end
                
    end
    
    tEnd=toc(tStart);
    
    if CtrlVar.InfoLevelNonLinIt>3 && iteration >= 2
        
        N=max([1,iteration-5]);
        
        [~,~,a1]=detrend_xt(log10(diffVector(N:iteration)),N:iteration);
        fprintf(CtrlVar.fidlog,' slope NR : %14.7g \n',a1);
        if CtrlVar.doplots==1
            figure ; semilogy(0:iteration,diffVector(1:iteration+1),'x-r') ;
            title('r^2 residual (NR uv diagnostic step)') ; xlabel('Iteration') ; ylabel('r^2 (Residual)')
        end
        
    end
    
    if isnan(r)
        fprintf(CtrlVar.fidlog,' SSTREAM2dNR returns NAN as residual!!! \n') ;
        warning('SSTREAM2NR:didnotconverge',' SSTREAM2dNR did not converge to a solution. Saving all variables in TestSaveNR.mat \n ')
        save TestSaveNR
    elseif r>CtrlVar.NLtol
        fprintf(CtrlVar.fidlog,' SSTREAM2dNR did NOT converge to given tolerance of %-g with r=%-g in %-i iterations and in %-g  sec \n',CtrlVar.NLtol,r,iteration,tEnd);
        RunInfo.converged=0;
        warning('SSTREAM2NR:didnotconverge',' SSTREAM2dNR did not converge to a solution. Saving all variables in TestSaveNR.mat \n ')
        save TestSaveNR
    else
        if CtrlVar.InfoLevelNonLinIt>0
            fprintf(CtrlVar.fidlog,' SSTREAM2dNR converged to given tolerance of %-g with r=%-g in %-i iterations and in %-g  sec \n',CtrlVar.NLtol,r,iteration,tEnd) ;
        end
    end
    
    
    if iteration > CtrlVar.NRitmax
        fprintf(CtrlVar.fidlog,'Maximum number of NR iterations %-i reached in uv loop with r=%-g \n',CtrlVar.NRitmax,r);
        warning('SSTREAM2dNR:MaxIterationReached','SSTREAM2NR exits because maximum number of iterations %-i reached with r=%-g \n',CtrlVar.NRitmax,r)
    end
    
    RunInfo.Iterations=iteration;  RunInfo.residual=r;
    
    if any(isnan(ub)) || any(isnan(vb))  ; save TestSaveNR  ;  error(' nan in ub vb ') ; end
    
    
   
end

