function [UserVar,RunInfo,F1,l1,BCs1]=SSTREAM_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)
    
    narginchk(8,8)
    nargoutchk(4,5)
    
    
    
    if CtrlVar.InfoLevelNonLinIt>=10  ; fprintf(CtrlVar.fidlog,' \n SSTREAM(uvh): Transient implicit with respect to u, v, and h  \n ') ; end
    
    %%
    % Fully implicit Newton-Raphson with regard to both u, v and h
    % advances the solution by dt
    %
    % h0, u0, and v0 are values at the start of the time step,
    % on input h1,u1,v1 are estimates for h, u, and v at the end of the time step
    % on exit  h1, u1, and v1 are calculated values for u,v and h at the end of the time step
    
    
    
    %%
    
    
    
    % I need to solve
    %
    % [Kxu Kxv Kxh Luv'  0  ] [du]        =  [ -Ru ] - Luv' luv
    % [Kyu Kyv Kyh          ] [dv]           [ -Rv ]
    % [Khu Khv Khh  0   Lh' ] [dh]           [ -Rh- Lh' lh ]
    % [  Luv        0    0  ] [duv]          [ cuv-Luv [u ;v ]
    % [ 0    Lh  0  0    0  ] [dlh]          [ ch-Lh h]
    %
    % All matrices are Nnodes x Nnodes, apart from:
    % Luv is #uv constraints x 2 Nnodes, i.e. Luv [u;v]= cuv
    % Lh  is # h contraints x Nnodes, i.e.    Lh h= ch
    %  or
    %
    % [K L'] [ duvh ]      =  [ -R- L' l ]
    % [L 0 ] [  dl  ]         [ cuvh-L [u;v;h]  ]
    %
    % where
    %
    % K= [Kxu Kxv Kxh]
    %    [Kyu Kyv Kyh]
    %    [Khu Khv Khh]
    %
    % and
    % L=[Luv 0]
    %   [0  Lh]
    % and uvh=[u;v;h], duvh=[du;dv; dh]  and l=[luv ; lh]
    % where L [u;v;h]=cuvh
    %
    
     
    rVector.gamma=zeros(CtrlVar.NRitmax+1,1)+NaN;
    rVector.ruv=zeros(CtrlVar.NRitmax+1,1)+NaN;
    rVector.rWork=zeros(CtrlVar.NRitmax+1,1)+NaN;
    rVector.rForce=zeros(CtrlVar.NRitmax+1,1)+NaN;
    rVector.D2=zeros(CtrlVar.NRitmax+1,1)+NaN;
    
    BackTrackSteps=0;
    
    
    if any(F0.h<0) ; warning('MATLAB:SSTREAM_TransientImplicit',' thickness negative ') ; end
    
    
    tStart=tic;
    
    
    [F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);  % make sure that if any extrapolation of fields or interpolation was performed,
    % that the geometrical fields are consistent with floation condition.
    % However, since the uvh formulation is with respect to h
    % alone, this will not affect the solution since this does not
    % change h.
    dub=F1.ub-F0.ub; dvb=F1.vb-F0.vb ; dh=F1.h-F0.h; 
    
    
    %%
    if CtrlVar.GuardAgainstWildExtrapolationInExplicit_uvh_Step
        N=3;
        
        speed1=sqrt(F1.ub.*F1.ub+F1.vb.*F1.vb);
        speed0=sqrt(F0.ub.*F0.ub+F0.vb.*F0.vb);
        Duv=(speed1-speed0)./(speed0+10*CtrlVar.SpeedZero);
        Dh=(F1.h-F0.h)./(F0.h+10*CtrlVar.ThickMin);
        
        
        
        %     figure ; histogram(Duv);
        %     figure ; histogram(Dh);
        
        
        Iuvh=((Duv-mean(Duv))> N*std(Duv)) | ((Dh-mean(Dh)) > N*std(Dh)) | abs(Duv)>0.1 | abs(Dh) > 0.1;
        
        
        fprintf(' Guarding agains wild extrapolation in uvh step.\n')
        fprintf(' Resetting %i forward explicit estimates out of %i to values at previous time step. \n',...
            numel(find(Iuvh)),numel(Iuvh))
        
        
        %     Iuvh=abs(dub-mean(dub))>N*std(dub) | abs(dvb-mean(dvb))>N*std(dvb) | abs(dh-mean(dh))>N*std(dh) ;
        %     fprintf(' Guarding agains wild extrapolation in uvh step.\n')
        %     fprintf(' Resetting %i forward explicit estimates out of %i to values at previous time step. \n',...
        %         numel(find(Iuvh)),numel(Iuvh))
        %
        
        
        
        F1.ub(Iuvh)=F0.ub(Iuvh);
        F1.vb(Iuvh)=F0.vb(Iuvh);
        F1.h(Iuvh)=F0.h(Iuvh);
    end
    
    
    
    
    
    
    
    
    %% assemble global Lagrange constraint matrix
    MLC=BCs2MLC(CtrlVar,MUA,BCs1);
    
    % Luv=MLC.ubvbL;
    % cuv=MLC.ubvbRhs;
    % Lh=MLC.hL;
    % ch=MLC.hRhs;
    
    if numel(l1.ubvb)~=numel(MLC.ubvbRhs) ; l1.ubvb=zeros(numel(MLC.ubvbRhs),1) ; end
    if numel(l1.h)~=numel(MLC.hRhs) ; l1.h=zeros(numel(MLC.hRhs),1) ; end
    nlubvb=numel(l1.ubvb) ;
    
    
    %[L,cuvh,luvh]=AssembleLuvh(Luv,Lh,cuv,ch,l1.ubvb,l1.h,MUA.Nnodes);
    [L,cuvh,luvh]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs1,l1);
    dl=luvh*0;
    
    if ~isempty(L)
        BCsRelativeError=norm(L*[F1.ub;F1.vb;F1.h]-cuvh)/norm(cuvh+1000*eps);
    else
        BCsRelativeError=0;
    end
    
    if BCsRelativeError>0.01
        
        fprintf('WARNING: At the beginning of the uvh iteration F1 is not a feasable point\n')
        % fprintf('         Although the uvh iteration can start at an infeasable point and still converge successfully, \n')
        
    end
    

    
    CtrlVar.uvhMatrixAssembly.ZeroFields=true;
    CtrlVar.uvhMatrixAssembly.Ronly=true;
    [UserVar,RunInfo,R0,~]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1);
    Fext0=R0;
    
    iteration=0 ;
    
    RunInfo.Forward.Converged=0;
    RunInfo.BackTrack.Converged=1 ;
    r=inf;  rWork=inf ; rForce=inf;
    gamma=1 ;
    
    
  
    
    while true
        
   
        
        % I define two types of exit critera:
        %  1) if residuals are below a given tolerance, exit and call it a success
        %  2) if step size is below a prescribed fraction of the Newton step in two
        %     consecutive iterations, the minimisation procedure is judged to have
        %     stagnated. There can be good or bad reasons for this: Possibly the
        %     minimisation procedure has converged on the miniumum and the norm of the
        %     gradient is small-enough (good) or it is too-large (bad).
        %
        
       
        
        
        if gamma > max(CtrlVar.uvhExitBackTrackingStepLength,CtrlVar.BacktrackingGammaMin)
            
            ResidualsCriteria=(rWork<CtrlVar.uvhDesiredWorkAndForceTolerances(1)  && rForce<CtrlVar.uvhDesiredWorkAndForceTolerances(2))...
                && (rWork<CtrlVar.uvhDesiredWorkOrForceTolerances(1)  || rForce<CtrlVar.uvhDesiredWorkOrForceTolerances(2))...
                && iteration >= CtrlVar.NRitmin;

            
        else
            
            ResidualsCriteria=(rWork<CtrlVar.uvhAcceptableWorkAndForceTolerances(1)  && rForce<CtrlVar.uvhAcceptableWorkAndForceTolerances(2))...
                && (rWork<CtrlVar.uvhAcceptableWorkOrForceTolerances(1)  || rForce<CtrlVar.uvhAcceptableWorkOrForceTolerances(2))...
                && iteration >= CtrlVar.NRitmin;
            
        end
        
        
        if ResidualsCriteria
            
            tEnd=toc(tStart);
            if CtrlVar.InfoLevelNonLinIt>=1
                fprintf(' SSTREAM(uvh) (time|dt)=(%g|%g): Converged with rForce=%-g and rWork=%-g in %-i iterations and in %-g  sec \n',...
                    CtrlVar.time,CtrlVar.dt,rForce,rWork,iteration,tEnd) ;
            end
            RunInfo.Forward.Converged=1;
            break
            
        end
        
        
        if iteration > CtrlVar.NRitmax
            
            if CtrlVar.InfoLevelNonLinIt>=1
                fprintf(' SSTREAM(uvh) (time|dt)=(%g|%g): Maximum number of non-linear iterations reached. uvh iteration did not converge! \n',CtrlVar.time,CtrlVar.dt)
                fprintf(' Exiting uvh iteration after %-i iterations with r=%-g \n',iteration,r)
            end
            
            RunInfo.Forward.Converged=0;
            break
        end
   
        if RunInfo.BackTrack.Converged==0 || gamma==0 
            if CtrlVar.InfoLevelNonLinIt>=1
                fprintf(' SSTREAM(uvh) (time|dt)=(%g|%g): Backtracting within non-linear iteration stagnated! \n Exiting non-lin iteration with r=%-g  after %-i iterations. \n',...
                    CtrlVar.time,CtrlVar.dt,r,iteration) ;
            end
            
            if CtrlVar.WriteRunInfoFile
                fprintf(RunInfo.File.fid,' SSTREAM(uvh) (time|dt)=(%g|%g): Backtracting within non-linear iteration stagnated! \n Exiting non-lin iteration with r=%-g   after %-i iterations. \n',...
                    CtrlVar.time,CtrlVar.dt,r,iteration) ;
            end
            
            RunInfo.Forward.Converged=0;
            break
        end
        
        iteration=iteration+1;
        
        
        %% Newton step
        % If I want to use the Newton Decrement (work) criterion I must calculate the Newton
        % step ahead of the cost function

        
        CtrlVar.uvhMatrixAssembly.ZeroFields=false; CtrlVar.uvhMatrixAssembly.Ronly=false;
        [UserVar,RunInfo,Ruvh,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1);
        
        if ~isempty(L)
            frhs=-Ruvh-L'*luvh;
            grhs=cuvh-L*[F1.ub;F1.vb;F1.h];
        else
            frhs=-Ruvh;
            grhs=[];
        end

        
        [duvh,dl]=solveKApe(K,L,frhs,grhs,[dub;dvb;dh],dl,CtrlVar);
        dub=duvh(1:MUA.Nnodes) ;  dvb=duvh(MUA.Nnodes+1:2*MUA.Nnodes); dh=duvh(2*MUA.Nnodes+1:end);
        
  
        
        % [r0,UserVar,RunInfo,rForce0,rWork0,D20]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,Fext0);
                                
        Func=@(gamma) CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,Fext0) ;
        gamma=0 ; [r0,UserVar,RunInfo,rForce0,rWork0,D20]=Func(gamma); 
        
        if iteration==1  % save the first r value for plotting, etc
            rVector.gamma(1)=gamma;
            rVector.ruv(1)=NaN;
            rVector.rWork(1)=rWork0;
            rVector.rForce(1)=rForce0 ;
            rVector.D2(1)=D20 ;
        end
        
        %% calculate  residuals at full Newton step, i.e. at gamma=1
  
           
        % r1=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,Fext0);
        gamma=1 ; [r1,UserVar,RunInfo,rForce1,rWork1,D21]=Func(gamma); 
        
        %% either accept full Newton step or do a line search
        
        % [UserVar,RunInfo,gamma,r]=FindBestGamma2DuvhBacktrack(UserVar,RunInfo,CtrlVar,MUA,F0,F1,dub,dvb,dh,dl,L,luvh,cuvh,r0,r1,Fext0);
        
        
        
        slope0=-2*r0 ; 
        [gamma,r,BackTrackInfo]=BackTracking(slope0,1,r0,r1,Func,CtrlVar);

        RunInfo.BackTrack=BackTrackInfo; 
        
        % If backtracking returns all values, then this call will not be needed.
        % [rTest,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,Fext0);
        [rTest,~,~,rForce,rWork,D2]=Func(gamma); 
        rVector.gamma(iteration+1)=gamma;
        rVector.ruv(iteration+1)=NaN;
        rVector.rWork(iteration+1)=rWork;
        rVector.rForce(iteration+1)=rForce ;
        rVector.D2(iteration+1)=D2 ;
        if ~isequal(r,rTest) 
            fprintf("r=%g \t rTest=%g \t \n",r,rTest)
            error('SSTREAM_TransientImplicit:expecting r and rTest to be equal')
        end
        
        
        %% If desired, plot residual along search direction
        if CtrlVar.InfoLevelNonLinIt>=10 && CtrlVar.doplots==1
            nnn=50;
            gammaTestVector=zeros(nnn,1) ; rForceTestvector=zeros(nnn,1);  rWorkTestvector=zeros(nnn,1); rD2Testvector=zeros(nnn,1);
            Upper=2.2;
            Lower=-1 ;
            if gamma>0.7*Upper ; Upper=2*gamma; end
            parfor I=1:nnn
                gammaTest=(Upper-Lower)*(I-1)/(nnn-1)+Lower
                [rTest,~,~,rForceTest,rWorkTest,D2Test]=Func(gammaTest);
                %[rTest,~,~,rForceTest,rWorkTest,D2Test]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gammaTest,Fext0);
                gammaTestVector(I)=gammaTest ; rForceTestvector(I)=rForceTest; rWorkTestvector(I)=rWorkTest;  rD2Testvector(I)=D2Test;
            end
            
            gammaZero=min(abs(gammaTestVector)) ;
            if gammaZero~=0
                [rTest,~,~,rForceTest,rWorkTest,D2Test]=Func(0);
                % [rTest,~,~,rForceTest,rWorkTest,D2Test]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,0,Fext0);
                gammaTestVector(nnn+1)=0 ; rForceTestvector(nnn+1)=rForceTest; rWorkTestvector(nnn+1)=rWorkTest;  rD2Testvector(nnn+1)=D2Test;
            end
            
            [gammaTestVector,ind]=unique(gammaTestVector) ; rForceTestvector=rForceTestvector(ind) ; rWorkTestvector=rWorkTestvector(ind) ; rD2Testvector=rD2Testvector(ind) ;
            [gammaTestVector,ind]=sort(gammaTestVector) ; rForceTestvector=rForceTestvector(ind) ; rWorkTestvector=rWorkTestvector(ind) ; rD2Testvector=rD2Testvector(ind) ;
            [temp,I0]=min(abs(gammaTestVector)) ;
            
            SlopeForce=-2*rForce0;
            SlopeWork=-2*rWork0;
            SlopeD2=-D20;
            CtrlVar.MinimisationQuantity=CtrlVar.uvhMinimisationQuantity;
            [ForceFig,WorkFig]=PlotCostFunctionsVersusGamma(CtrlVar,RunInfo,gamma,r,iteration,"-uvh-",...
                gammaTestVector,rForceTestvector,rWorkTestvector,rD2Testvector,...
                SlopeForce,SlopeWork,SlopeD2,rForce,rWork,D2);

        end
        
        
   
        %% update variables
        
        
        F1.ub=F1.ub+gamma*dub;
        F1.vb=F1.vb+gamma*dvb;
        F1.h=F1.h+gamma*dh;
        luvh=luvh+gamma*dl;
        
        l1.ubvb=luvh(1:nlubvb) ;  l1.h=luvh(nlubvb+1:end);
        
        
        
        temp=CtrlVar.ResetThicknessToMinThickness;
        if ~CtrlVar.ResetThicknessInNonLinLoop
            CtrlVar.ResetThicknessToMinThickness=0;
        end
        
        % make sure to update s and b as well!
        [F1.b,F1.s]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);  
        CtrlVar.ResetThicknessToMinThickness=temp;
 
        
        if ~isempty(L)
            BCsError=norm(L*[F1.ub;F1.vb;F1.h]-cuvh);
        else 
            BCsError=0;
        end
        
        % Variables have been updated, if I have MassBalanceGeometryFeedback>0 I must
        % update the surface mass balance within this non-linear loop. Actually I here
        % only need to consider option 1 because if options 2 or 3 are used the
        % mass-blance is updated anyhow witin the assmebly loop.
        if CtrlVar.MassBalanceGeometryFeedback>0
            
            rdamp=CtrlVar.MassBalanceGeometryFeedbackDamping;
            if rdamp~=0
                as1Old=F1.as ; ab1Old=F1.ab;
            end
            CtrlVar.time=CtrlVar.time+CtrlVar.dt;
            [UserVar,F1]=GetMassBalance(UserVar,CtrlVar,MUA,F1);
            CtrlVar.time=CtrlVar.time-CtrlVar.dt;
            
            if rdamp~=0
                % If Hessian inaccurate, or too non-linear, then dampen these changes might be a
                % good idea.
                F1.as=(1-rdamp)*F1.as+rdamp*as1Old;
                F1.ab=(1-rdamp)*F1.ab+rdamp*ab1Old;
            end
        end

        
        if CtrlVar.InfoLevelNonLinIt>=100  && CtrlVar.doplots==1
            PlotForceResidualVectors('uvh',Ruvh,L,luvh,MUA.coordinates,CtrlVar) ; axis equal tight
        end
     
        if CtrlVar.InfoLevelNonLinIt>=1
            fprintf(...
                'NR-SSTREAM(uvh):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , rForce=%-14.7g , rWork=%-14.7g , BCsError=%-g  \n ',...
                iteration,RunInfo.BackTrack.iarm,gamma,r/r0,r0,r,rForce,rWork,BCsError);
            
        end
        
        
        
        if CtrlVar.WriteRunInfoFile
            
            fprintf(RunInfo.File.fid,...
                'NR-SSTREAM(uvh):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , rForce=%-14.7g , rWork=%-14.7g , BCsError=%-g   \n ',...
                iteration,RunInfo.BackTrack.iarm,gamma,r/r0,r0,r,rForce,rWork,BCsError);
            
        end
        
        BackTrackSteps=BackTrackSteps+RunInfo.BackTrack.iarm ; 
        
        
    end
    
    %% return calculated values at the end of the time step
    %F1.ub=ub ; F1.vb=vb ; F1.h=h; l1.ubvb=luv1  ; l1.h=lh;
    
    % I got out of the while loop if either if the solver converged, or
    % backtrack stagnated.
    %RunInfo.Forward.Converged=1;
    if RunInfo.BackTrack.Converged==0
        RunInfo.Forward.Converged=0;
    end
    
    %% print/plot some info
    
    if CtrlVar.InfoLevelNonLinIt>=10 && iteration >= 2 && CtrlVar.doplots==1
        
        FindOrCreateFigure("NR-uvh r");
        yyaxis left
        semilogy(0:iteration,rVector.rForce(1:iteration+1),'x-') ;
        ylabel('rForce^2')
        yyaxis right
        semilogy(0:iteration,rVector.rWork(1:iteration+1),'o-') ;
        ylabel('rWork^2')
        
        title('Force and Work residuals (NR uvh diagnostic step)') ; xlabel('Iteration') ;
        
        
    end
    
    if ~isempty(L)
        if BCsError>10*eps
            fprintf(CtrlVar.fidlog,'Norm of BCs residuals is %14.7g  \n ',BCsError);
        end
    end
    
    
    tEnd=toc(tStart);
    
    
    
    if iteration > CtrlVar.NRitmax
        fprintf(CtrlVar.fidlog,'Warning: maximum number of NRuvh iterations %-i reached \n',CtrlVar.NRitmax);
        warning('SSTREAM2dNR:MaxIterationReached','SSTREAM2NR exits because maximum number of iterations %-i reached \n',CtrlVar.NRitmax)
        filename='Dumpfile_SSTREAM_TransientImplicit.mat';
        fprintf('Saving all data in a dumpfile %s \n',filename)
        save(filename)
    end
    
    RunInfo.Forward.Iterations=iteration; % May try to get rid of and use the uvhIterations vector instead
    RunInfo.Forward.uvhIterations(CtrlVar.CurrentRunStepNumber)=iteration ; 
    RunInfo.Forward.uvhResidual(CtrlVar.CurrentRunStepNumber)=r;
    RunInfo.Forward.uvhBackTrackSteps(CtrlVar.CurrentRunStepNumber)=BackTrackSteps ; 
    
    if CtrlVar.WriteRunInfoFile
        
        fprintf(RunInfo.File.fid,' --->  SSTREAM(uvh/%s) \t time=%15.5f \t dt=%-g \t r=%-g \t #it=% i \t CPUsec=%-g \n',...
            CtrlVar.uvhImplicitTimeSteppingMethod,CtrlVar.time,CtrlVar.dt,RunInfo.Forward.Residual,RunInfo.Forward.Iterations,tEnd) ;
        
    end
    
    
end


