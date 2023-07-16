function  [UserVar,F,l,Kuv,Ruv,RunInfo,L]=SSTREAM2dNR2(UserVar,CtrlVar,MUA,BCs,F,l,RunInfo)
    
    
    % Solves SSA/SSTREAM for u and v


    nargoutchk(7,7)
    narginchk(7,7)
    
    
    tStart=tic;
    RunInfo.Forward.uvConverged=1; 
 
    
    if isempty(CtrlVar.CurrentRunStepNumber) || CtrlVar.CurrentRunStepNumber==0 
        CtrlVar.CurrentRunStepNumber=1;
    end

   % RunInfo.Forward.uvIterations(CtrlVar.CurrentRunStepNumber)=NaN;  
   % RunInfo.Forward.Residual=NaN; BackTrackInfo.iarm=NaN;
    
    Kuv=[] ; Ruv=[]; 
   
    
    % MLC=BCs2MLC(CtrlVar,MUA,BCs);
    % L=MLC.ubvbL;
    % cuv=MLC.ubvbRhs;
    
    [L,cuv]=AssembleLuvSSTREAM(CtrlVar,MUA,BCs) ;
    
    
    if isempty(cuv)
        l.ubvb=[];
    elseif numel(l.ubvb)~=numel(cuv)
        l.ubvb=zeros(numel(cuv),1) ;
    end
    
    
    if any(isnan(F.C))
        save TestSave ;
        warning('SSTREAM2NR:CisNaN',' nan in C. Returning with NaN in solution.\n ') ;
        F.ub=F.ub+NaN;
        F.vb=F.vb+NaN;
        return
    end
    
    if any(isnan(F.AGlen))
        save TestSave
        warning('SSTREAM2NR:CisNaN',' nan in A. Returning with NaN in solution.\n ') ;
        F.ub=F.ub+NaN;
        F.vb=F.vb+NaN;
        return
    end
    
    
    
    if any(isnan(F.S)) ; save TestSave ; error( ' S nan ') ; end
    if any(isnan(F.h)) ; save TestSave  ; error( ' h nan ') ; end
    if any(isnan(F.ub)) ; save TestSave ; error( ' ub nan ') ; end
    if any(isnan(F.vb)) ; save TestSave ; error( ' vb nan ') ; end
    if any(isnan(l.ubvb)) ; save TestSave ; error( ' ubvbLambda nan ') ; end
    if any(isnan(F.rho)) ; save TestSave  ; error( ' rho nan ') ; end
    if any(F.h<0) ; warning('MATLAB:SSTREAM2dNR:hnegative',' thickness negative ') ; end
    
    
    
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
    
    % I write the system as
    % [Kuv  Luv^T ]  [duv]  =  [ -R(uv) - Luv^T l]
    % [Luv   0    ]  [dl]      [cuv-Luv uv]
    %
    
    
    %% Make sure iterate is feasable
    F.ub(BCs.ubFixedNode)=BCs.ubFixedValue; F.vb(BCs.vbFixedNode)=BCs.vbFixedValue;
    %%
    
    
    dub=zeros(MUA.Nnodes,1) ; dvb=zeros(MUA.Nnodes,1) ; dl=zeros(numel(l.ubvb),1);
    
    
    % diffVector=zeros(CtrlVar.NRitmax+1,1)+NaN;
    rVector.gamma=zeros(CtrlVar.NRitmax+1,1)+NaN;
    rVector.rDisp=zeros(CtrlVar.NRitmax+1,1)+NaN;
    rVector.rWork=zeros(CtrlVar.NRitmax+1,1)+NaN;
    rVector.rForce=zeros(CtrlVar.NRitmax+1,1)+NaN;
    
    Kuv=[] ; 
    
     
    
 

    fext0=KRTFgeneralBCs(CtrlVar,MUA,F,true); % RHS with velocities set to zero, i.e. only external forces
    
    %% New normalisation idea, 10 April 2023
    % set (ub,vb) to zero, except where BCs imply otherwise, ie make the iterate feasable 
    % then calculate the const function for this value and use as normalisation
    % gamma=0; fext0=1; 
    % ubStart=F.ub; vbStart=F.vb; 
    % F.ub=zeros(MUA.Nnodes,1); F.vb=zeros(MUA.Nnodes,1); 
    % F.ub(BCs.ubFixedNode)=BCs.ubFixedValue; F.vb(BCs.vbFixedNode)=BCs.vbFixedValue; % Make sure iterate is feasable
    % fNOrm=CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,dub,dvb,dl) ;
    % fext0=sqrt(fNOrm); 
    % F.ub=ubStart ; F.vb=vbStart ; 
    %%

    
    
    Ruv=KRTFgeneralBCs(CtrlVar,MUA,F);     % RHS with calculated velocities, i.e. difference between external and internal forces
    
    RunInfo.CPU.Solution.uv=0;

 

    % The initial estimate must be based on residuals as both displacements and work
    % requires solving the Newton system.
    % gamma=0 ; [r,UserVar,RunInfo,rForce,rWork] = CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,dub,dvb,dl) ;
    Func=@(gamma) CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,dub,dvb,dl) ;
    gamma=0 ; [r,UserVar,RunInfo,rForce]=Func(gamma);
    rWork=0 ; % I set rWork to zero to disable it as initial criterion
    gamma=1;
    
    iteration=0;  ResidualReduction=1e10; RunInfo.CPU.solution.uv=0 ; RunInfo.CPU.Assembly.uv=0;

    BackTrackInfo.iarm=nan;
    while true

        
        ResidualsCriteria=uvResidualsCriteria(CtrlVar,rForce,rWork,iteration,gamma) ; 

        

        if ResidualsCriteria
            
            tEnd=toc(tStart);
            if CtrlVar.InfoLevelNonLinIt>=1
                fprintf(' SSTREAM(uv) (time|dt)=(%g|%g): Converged with rForce=%-g and rWork=%-g in %-i iterations and in %-g  sec \n',...
                    CtrlVar.time,CtrlVar.dt,rForce,rWork,iteration,tEnd) ;
            end
            RunInfo.Forward.uvConverged=1;
            break
            
        end

        
        if iteration > CtrlVar.NRitmax
            
            if CtrlVar.InfoLevelNonLinIt>=1
                fprintf(' SSTREAM(uv) (time|dt)=(%g|%g): Maximum number of non-linear iterations reached. uv iteration did not converge! \n',CtrlVar.time,CtrlVar.dt)
                fprintf(' Exiting uv iteration after %-i iterations with r=%-g \n',iteration,r)
            end
            
            if CtrlVar.WriteRunInfoFile
                fprintf(RunInfo.File.fid,' SSTREAM(uv) (time|dt)=(%g|%g): Maximum number of non-linear iterations reached. uv iteration did not converge! \n',CtrlVar.time,CtrlVar.dt);
                fprintf(RunInfo.File.fid,' Exiting uv iteration after %-i iterations with r=%-g \n',iteration,r);
            end
            
            RunInfo.Forward.uvConverged=0; 
            break
        end
        
        
        iteration=iteration+1;
        
        %% Newton step
        % If I want to use the Newton Decrement (work) criterion I must calculate the Newton
        % step ahead of the cost function
        if rem(iteration-1,CtrlVar.ModifiedNRuvIntervalCriterion)==0  || ResidualReduction> CtrlVar.ModifiedNRuvReductionCriterion
            
            tAssembly=tic;
            [Ruv,Kuv]=KRTFgeneralBCs(CtrlVar,MUA,F);
            RunInfo.CPU.Assembly.uv=toc(tAssembly)+RunInfo.CPU.Assembly.uv;
            NRincomplete=0;
        else
            tAssembly=tic;
            Ruv=KRTFgeneralBCs(CtrlVar,MUA,F);
            RunInfo.CPU.Assembly.uv=toc(tAssembly)+RunInfo.CPU.Assembly.uv;
            NRincomplete=1;
        end

        
        if CtrlVar.TestAdjointFiniteDifferenceType=="complex step differentiation"
            CtrlVar.TestForRealValues=false;
        end
        
        if CtrlVar.TestForRealValues
            if ~isreal(Kuv) ; save TestSave Kuv ; error('SSTREAM2dNR: K not real') ;  end
            if ~isreal(L) ; save TestSave L ; error('SSTREAM2dNR: L not real') ;  end
        end
        
        if any(isnan(Kuv)) ; save TestSave Kuv ; error('SSTREAM2dNR: K nan') ;  end
        if any(isnan(L)) ; save TestSave L ; error('SSTREAM2dNR: L nan') ;  end
        
        CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical=issymmetric(Kuv) ;

        
        % [Kuv  Luv' ]  [duv]  =  [ -R(uv) - Luv' l]
        % [Luv   0    ]  [dl]      [ cuv-Luv uv      ]
        
        if ~isempty(L)
            frhs=-Ruv-L'*l.ubvb;
            grhs=cuv-L*[F.ub;F.vb];
        else
            frhs=-Ruv;
            grhs=[];
        end
   
        tSolution=tic;
        
        if CtrlVar.Solver.isUpperLeftBlockMatrixSymmetrical
            [sol,dl]=solveKApeSymmetric(Kuv,L,frhs,grhs,[dub;dvb],dl,CtrlVar);
        else
            [sol,dl]=solveKApe(Kuv,L,frhs,grhs,[dub;dvb],dl,CtrlVar);
        end
        
        RunInfo.CPU.Solution.uv=toc(tSolution)+RunInfo.CPU.Solution.uv;
        
        
        if CtrlVar.TestForRealValues
            dub=real(sol(1:MUA.Nnodes)) ; dvb=real(sol(MUA.Nnodes+1:2*MUA.Nnodes));
        else
            dub=sol(1:MUA.Nnodes) ; dvb=sol(MUA.Nnodes+1:2*MUA.Nnodes);
        end

        %% Residuals , at gamma=0;
        Func=@(gamma) CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,dub,dvb,dl) ;
        gamma=0 ; [r0,UserVar,RunInfo,rForce0,rWork0,D20]=Func(gamma);
        
        if iteration==1
%           special case to check if initial state was already within tolerance 
            ResidualsCriteria=uvResidualsCriteria(CtrlVar,rForce0,rWork0,iteration,1) ;  % here I set gamma as an input to one, i.e. as if I had take a full step
            
            % save the first r value for plotting, etc
            rVector.gamma(1)=gamma;
            rVector.rDisp(1)=NaN;
            rVector.rWork(1)=rWork0;
            rVector.rForce(1)=rForce0 ;
            if ResidualsCriteria
                iteration=iteration-1 ; % reset to zero as the iteration was never performed 
                tEnd=toc(tStart);
                if CtrlVar.InfoLevelNonLinIt>=1
                    fprintf(' SSTREAM(uv) (time|dt)=(%g|%g): Converged with rForce=%-g and rWork=%-g in %-i iterations and in %-g  sec \n',...
                        CtrlVar.time,CtrlVar.dt,rForce0,rWork0,iteration,tEnd) ;
                end
                RunInfo.Forward.uvConverged=1;

                break
            end
        end



        %% calculate  residuals at full Newton step, i.e. at gamma=1
        gamma=1 ; [r1,UserVar,RunInfo,rForce1,rWork1,D21]=Func(gamma);

        if r1/r0 < CtrlVar.NewtonAcceptRatio

            r=r1 ; rForce=rForce1 ; rWork=rWork1 ; D2=D21 ;
            du=dub ; dv=dvb ;
            BackTrackInfo.Infovector=[0 r0 ; 1 r1] ;
            BackTrackInfo.Converged=1; BackTrackInfo.iarm=0; 

        else


            func=@(gamma,Du,Dv,Dl) CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,Du,Dv,Dl) ;
            dh=[] ; dJdh=[] ;
            dJdu=frhs(1:MUA.Nnodes);
            dJdv=frhs(MUA.Nnodes+1:2*MUA.Nnodes);
            dJdl=grhs ;
            Normalisation=fext0'*fext0+1000*eps;
            % CtrlVar.InfoLevelBackTrack=1000;  CtrlVar.InfoLevelNonLinIt=10 ;
            [gamma,r,du,dv,dh,dl,BackTrackInfo,rForce,rWork,D2] = rLineminUa(CtrlVar,UserVar,func,r0,r1,Kuv,L,dub,dvb,dh,dl,dJdu,dJdv,dJdh,dJdl,Normalisation,MUA.M) ;


        end


        RunInfo.BackTrack=BackTrackInfo;

        rVector.gamma(iteration+1)=gamma;
        rVector.rDisp(iteration+1)=NaN;
        rVector.rWork(iteration+1)=rWork;
        rVector.rForce(iteration+1)=rForce ;



        if BackTrackInfo.Converged==0
            fprintf(CtrlVar.fidlog,' SSTREAM2dNR backtracking step did not converge \n ') ;
            warning('SSTREAM2NR:didnotconverge',' SSTREAM2dNR backtracking step did not converge \n ')
            fprintf(CtrlVar.fidlog,' saving variables in SSTREAM2dNRDump \n ') ;
            save SSTREAM2dNRDump
            RunInfo.Forward.uvConverged=0; 
            break
        end


    
        %%

        %% If requested, plot residual as function of steplength
        if CtrlVar.InfoLevelNonLinIt>=10 && CtrlVar.doplots==1
            nnn=50;
            gammaTestVector=zeros(nnn,1) ; rForceTestvector=zeros(nnn,1); rWorkTestvector=zeros(nnn,1); rD2Testvector=zeros(nnn,1);
            
            Up=2;
            if gamma>0.7*Up ; Up=2*gamma; end
            parfor I=1:nnn
                gammaTest=Up*(I-1)/(nnn-1)+gamma/1000;
                [rTest,~,~,rForceTest,rWorkTest,D2Test]=Func(gammaTest);
                gammaTestVector(I)=gammaTest ; rForceTestvector(I)=rForceTest; rWorkTestvector(I)=rWorkTest; rD2Testvector(I)=D2Test; 
            end
            
            [gammaTestVector,ind]=unique(gammaTestVector) ; rForceTestvector=rForceTestvector(ind) ; rWorkTestvector=rWorkTestvector(ind) ;  rD2Testvector=rD2Testvector(ind) ;
            [gammaTestVector,ind]=sort(gammaTestVector) ; rForceTestvector=rForceTestvector(ind) ; rWorkTestvector=rWorkTestvector(ind) ; rD2Testvector=rD2Testvector(ind) ;
            
            SlopeForce=-2*rForce0;   
            SlopeWork=-2*rWork0;
            SlopeD2=-D20;
            CtrlVar.MinimisationQuantity=CtrlVar.uvMinimisationQuantity;
            PlotCostFunctionsVersusGamma(CtrlVar,RunInfo,gamma,r,iteration,"-uv-",...
                gammaTestVector,rForceTestvector,rWorkTestvector,rD2Testvector,...
                SlopeForce,SlopeWork,SlopeD2,rForce,rWork,D2);
            
        end

        % Need to update all primary (u,v,l) and dependent variables.
        % Here I have no dependent variables
        % F.ub=F.ub+gamma*dub ;
        % F.vb=F.vb+gamma*dvb;
        % l.ubvb=l.ubvb+gamma*dl;

        F.ub=F.ub+du ;
        F.vb=F.vb+dv;
        l.ubvb=l.ubvb+dl;
        
        ResidualReduction=r/r0;
        
        if CtrlVar.InfoLevelNonLinIt>100  && CtrlVar.doplots==1
            %PlotForceResidualVectors2('uv',Ruv,L,l.ubvb,MUA.coordinates,CtrlVar) ; axis equal tight
            PlotForceResidualVectors2(CtrlVar,MUA,F,'uv',Ruv,L,l.ubvb,iteration);
        end
        
        if NRincomplete
            stri='i';
        else
            stri=[];
        end
        
        if CtrlVar.InfoLevelNonLinIt>=1
            
            fprintf(CtrlVar.fidlog,'%sNR-SSTREAM(uv):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , rForce=%-14.7g , rWork=%-14.7g , fSlope0=%-14.7g \n ',...
                stri,iteration,RunInfo.BackTrack.iarm,gamma,r/r0,r0,r,rForce,rWork,2*rForce0);
        end
        
        if CtrlVar.WriteRunInfoFile
            
            fprintf(RunInfo.File.fid,'%sNR-SSTREAM(uv):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , rWork=%-14.7g , Assembly=%f sec. Solution=%f sec.\n ',...
                stri,iteration,RunInfo.BackTrack.iarm,gamma,r/r0,r0,r,rWork,RunInfo.CPU.Assembly.uv,RunInfo.CPU.Solution.uv);
        end
        
    end
    
    if CtrlVar.InfoLevelNonLinIt>=10 && iteration >= 2 && CtrlVar.doplots==1
  
            
            figruv=FindOrCreateFigure("NR-uv r"); clf(figruv) ;
            yyaxis left
            semilogy(0:iteration,rVector.rForce(1:iteration+1),'x-') ;
            ylabel('rResiduals^2')
            yyaxis right
            semilogy(0:iteration,rVector.rWork(1:iteration+1),'o-') ;
            ylabel('rWork^2')
            
            title('Force and Work residuals (NR uv diagnostic step)') ; xlabel('Iteration') ;
            
  
    end
    
    if isnan(r)
        fprintf(CtrlVar.fidlog,' SSTREAM2dNR returns NAN as residual!!! \n') ;
        warning('uvSSTREAM:didnotconverge',' SSTREAM2dNR did not converge to a solution. Saving all variables in TestSaveNR.mat \n ')
        save TestSaveNR
        
    end


    if numel(RunInfo.Forward.uvIterations) < CtrlVar.CurrentRunStepNumber
        RunInfo.Forward.uvIterations=[RunInfo.Forward.uvIterations;RunInfo.Forward.uvIterations+NaN];
        RunInfo.Forward.uvResidual=[RunInfo.Forward.uvResidual;RunInfo.Forward.uvResidual+NaN];
        RunInfo.Forward.uvBackTrackSteps=[RunInfo.Forward.uvBackTrackSteps;RunInfo.Forward.uvBackTrackSteps+NaN];
    end
    
    RunInfo.Forward.uvIterations(CtrlVar.CurrentRunStepNumber)=iteration;  
    RunInfo.Forward.uvResidual(CtrlVar.CurrentRunStepNumber)=r;
    RunInfo.Forward.uvBackTrackSteps(CtrlVar.CurrentRunStepNumber)=BackTrackInfo.iarm ; 


    
    if any(isnan(F.ub)) || any(isnan(F.vb))  ; save TestSaveNR  ;  error(' nan in ub vb ') ; end
    
    
    
end

