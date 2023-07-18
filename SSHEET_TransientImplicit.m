function [UserVar,RunInfo,F1,l1,BCs]=SSHEET_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs)

%  s and h are the initial estimates for s1 and h1
%  these are then updated, once convergent, I set s1=s and h1=h


nargoutchk(5,5)
narginchk(8,8)



[F0.b,F0.s,F0.h,F0.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F0.h,F0.S,F0.B,F0.rho,F0.rhow);
[F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);



a0=F0.as+F0.ab ; a1=F1.as+F1.ab;
if CtrlVar.InfoLevelNonLinIt>=10  ; fprintf(CtrlVar.fidlog,' \n uvh2DSSHEET \n ') ; end

%%
% SSHEET implicit Newton-Raphson with regard to  h
% advances the solution by dt
%
% The natural boundary conditions is no-flux.
%
% h0 is value at the start of the time step,
% on input h1 is an estimate for h at the end of the time step
% on exit  h1, u1, and v1 are calculated values for u,v and h at the end of the time step


% I need to solve
% [Khh Lh'] [ dh ]      =  [ -Rh-Lh' lambdah ]
% [Lh   0 ] [dlambdah]     [  Lhrhs-Lh h]


rVector.gamma=zeros(CtrlVar.NRitmax+1,1)+NaN;
rVector.ruv=zeros(CtrlVar.NRitmax+1,1)+NaN;
rVector.rWork=zeros(CtrlVar.NRitmax+1,1)+NaN;
rVector.rForce=zeros(CtrlVar.NRitmax+1,1)+NaN;
rVector.D2=zeros(CtrlVar.NRitmax+1,1)+NaN;



if any(F0.h<0) ; warning('MATLAB:uvh2DSSHEET',' thickness negative ') ; end

tStart=tic;

MLC=BCs2MLC(CtrlVar,MUA,BCs);
Lh=MLC.hL ; ch=MLC.hRhs ;

if numel(l1.h)~=numel(ch) ; l1.h=zeros(numel(ch),1) ; end

dlambdah=l1.h*0;
dh=F0.h*0;
iteration=0 ; 

gamma=1;  rWork=inf ; rForce=inf ; 

while true
    
       if gamma > max(CtrlVar.hExitBackTrackingStepLength,CtrlVar.BacktrackingGammaMin)
            
            ResidualsCriteria=(rWork<CtrlVar.hDesiredWorkAndForceTolerances(1)  && rForce<CtrlVar.hDesiredWorkAndForceTolerances(2))...
                && (rWork<CtrlVar.hDesiredWorkOrForceTolerances(1)  || rForce<CtrlVar.hDesiredWorkOrForceTolerances(2))...
                && iteration >= CtrlVar.NRitmin;

        else
            
            ResidualsCriteria=(rWork<CtrlVar.hAcceptableWorkAndForceTolerances(1)  && rForce<CtrlVar.hAcceptableWorkAndForceTolerances(2))...
                && (rWork<CtrlVar.hAcceptableWorkOrForceTolerances(1)  || rForce<CtrlVar.hAcceptableWorkOrForceTolerances(2))...
                && iteration >= CtrlVar.NRitmin;
            
        end
        
        if ResidualsCriteria
            
            tEnd=toc(tStart);
            if CtrlVar.InfoLevelNonLinIt>=1
                fprintf(' SSHEET(h) (time|dt)=(%g|%g): Converged with rForce=%-g and rWork=%-g in %-i iterations and in %-g  sec \n',...
                    CtrlVar.time,CtrlVar.dt,rForce,rWork,iteration,tEnd) ;
            end
            RunInfo.Forward.uvhConverged=1;  RunInfo.Forward.hConverged=1;
            
            break
            
        end

        
        if iteration > CtrlVar.NRitmax
            
            if CtrlVar.InfoLevelNonLinIt>=1
                fprintf(' SSHEET(h) (time|dt)=(%g|%g): Maximum number of non-linear iterations reached. uvh iteration did not converge! \n',CtrlVar.time,CtrlVar.dt)
                fprintf(' Exiting h iteration after %-i iterations with r=%-g \n',iteration,r)
            end
            
            if CtrlVar.WriteRunInfoFile
                fprintf(RunInfo.File.fid,' SSHEET(h) (time|dt)=(%g|%g): Maximum number of non-linear iterations reached. uvh iteration did not converge! \n',CtrlVar.time,CtrlVar.dt);
                fprintf(RunInfo.File.fid,' Exiting h iteration after %-i iterations with r=%-g \n',iteration,r);
            end
            
            RunInfo.Forward.uvhConverged=0;  RunInfo.Forward.hConverged=0;
            
            break
        end

    
    iteration=iteration+1;
    
    
    if ~isequal(F1.dt,CtrlVar.dt) ; error("Ua:SSHEET_TransientImplict","F.dt not equal to CtrlVar.dt") ; end 
    [R,K,FF,T]=MatrixAssemblySSHEETtransient2HD(CtrlVar,MUA,F1.AGlen,F1.n,F1.C,F1.m,F1.rho,F1.g,F0.h,F0.b,F1.h,F1.b,a0,a1,F1.dt);
    
    if iteration==1
        FF0=FF ;  % There is a potential issue here which is that F0 is zero if the accumulation
                % and grad q is everywhere zero. However, if this happens dh/dt will also automatically be
                % zero. 
                % F0 is used as a normalisation factor when calculating the residual,
                % do not change this normalisation factor in the course of the iteration
    end

    
 
    
    %% solve the linear system
    if ~isempty(Lh)
        frhs=-R-Lh'*l1.h;
        grhs=ch-Lh*F1.h;
    else
        frhs=-R;
        grhs=[];
    end
    
    
    [dh,dlambdah]=solveKApe(K,Lh,frhs,grhs,dh,dlambdah,CtrlVar);
    
    
    if any(isnan(dh)) ; save TestSave  ;
        fprintf(CtrlVar.fidlog,'error: NaN in solution of implicit system \n') ;
        error(' NaN in solution of implicit h system ' ) ;
    end
    
    
    %% calculate  residuals at beginnin and for full Newton step
    
    Func=@(gamma) CalcCostFunctionSSHEET(UserVar,RunInfo,CtrlVar,gamma,dh,MUA,F1.AGlen,F1.n,F1.C,F1.m,F1.rho,F1.g,F0.h,F0.b,F1.h,F1.b,a0,a1,F1.dt,Lh,l1.h,dlambdah,FF0,ch);
    
    gamma=0; [r0,~,~,rForce0,rWork0,D20]=Func(gamma);
    gamma=1; [r1,~,~,rForce1,rWork1,D21]=Func(gamma);
    
    if iteration==1  % save the first r value for plotting, etc
        rVector.gamma(1)=gamma;
        rVector.ruv(1)=NaN;
        rVector.rWork(1)=rWork0;
        rVector.rForce(1)=rForce0 ;
        rVector.D2(1)=D20 ;
    end
    
    %% either accept full Newton step or do a line search
    
    Slope0=-2*r0 ;  % using the inner product def
    [gamma,r,BackTrackInfo]=BackTracking(Slope0,1,r0,r1,Func,CtrlVar);
    
    RunInfo.BackTrack=BackTrackInfo;
    
    
    
    [rTest,~,~,rForce,rWork,D2]=Func(gamma);
    rVector.gamma(iteration+1)=gamma;
    rVector.ruv(iteration+1)=NaN;
    rVector.rWork(iteration+1)=rWork;
    rVector.rForce(iteration+1)=rForce ;
    rVector.D2(iteration+1)=D2 ;
    if ~isequal(r,rTest)
        fprintf("r=%g \t rTest=%g \t \n",r,rTest)
        error('SSHEET_TransientImplicit:expecting r and rTest to be equal')
    end
    
    %% If desired, plot residual along search direction
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
        CtrlVar.MinimisationQuantity=CtrlVar.hMinimisationQuantity;
        [ForceFig,WorkFig]=PlotCostFunctionsVersusGamma(CtrlVar,RunInfo,gamma,r,iteration,"-h-",...
            gammaTestVector,rForceTestvector,rWorkTestvector,rD2Testvector,...
            SlopeForce,SlopeWork,SlopeD2,rForce,rWork,D2);
        
    end
    
    
   
    
    %% update variables
    F1.h=F1.h+gamma*dh ; l1.h=l1.h+gamma*dlambdah;
    
    temp=CtrlVar.ResetThicknessToMinThickness;
    if ~CtrlVar.ResetThicknessInNonLinLoop
        CtrlVar.ResetThicknessToMinThickness=0;
    end
    [F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);
    CtrlVar.ResetThicknessToMinThickness=temp;
    %[b1,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);



    if~isempty(Lh)
        BCsNorm=norm(ch-Lh*F1.h);
    else
        BCsNorm=0;
    end


    
    %% plot and print info
    if CtrlVar.InfoLevelNonLinIt>=100  && CtrlVar.doplots==1
        PlotForceResidualVectors('h-only',R,Lh,l1.h,MUA.coordinates,CtrlVar) ; axis equal tight
    end
    
    
    if CtrlVar.InfoLevelNonLinIt>=1
        fprintf(CtrlVar.fidlog,'NR-SSHEET(h):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , rForce=%-14.7g , rWork=%-14.7g , BCsNorm=%-14.7g\n ',...
            iteration,BackTrackInfo.iarm,gamma,r/r0,r0,r,rForce,rWork,BCsNorm);
    end
    
    
end

%% return calculated values at the end of the time step
% F1.h=h ;  lambdah1=l1.h;
[~,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);

[F1.ud,F1.vd,F1.ub,F1.vb]=uvSSHEET(CtrlVar,MUA,BCs,F1.AGlen,F1.n,F1.C,F1.m,F1.rho,F1.g,F1.s,F1.h);

tEnd=toc(tStart);

RunInfo.Forward.uvhIterations(CtrlVar.CurrentRunStepNumber)=iteration ;

if RunInfo.BackTrack.Converged==0
    RunInfo.Forward.uvhConverged=0;  RunInfo.Forward.hConverged=0;
end


%% print/plot some info
if CtrlVar.InfoLevelNonLinIt>=1
    fprintf(CtrlVar.fidlog,' NR-SSHEET(h) converged to given tolerance with r=%-g in %-i iterations and in %-g  sec \n',r,iteration,tEnd) ;
end

if CtrlVar.InfoLevelNonLinIt>=10 && iteration >= 2 && CtrlVar.doplots==1
    
    FindOrCreateFigure("NR-h r");
    yyaxis left
    semilogy(0:iteration,rVector.rForce(1:iteration+1),'x-') ;
    ylabel('rForce^2')
    yyaxis right
    semilogy(0:iteration,rVector.rWork(1:iteration+1),'o-') ;
    ylabel('rWork^2')
    
    title('Force and Work residuals (NR h SHEET diagnostic step)') ; xlabel('Iteration') ;
    
    
end


if iteration > CtrlVar.NRitmax
    fprintf(CtrlVar.fidlog,'Warning: maximum number of NRh iterations %-i reached \n',CtrlVar.NRitmax);
    warning('uvh2DSSHEET:MaxIterationReached','uvh2DSSHEET exits because maximum number of iterations %-i reached \n',CtrlVar.NRitmax)
    RunInfo.Forward.hConverged=0;
end

RunInfo.Forward.hIterations(CtrlVar.CurrentRunStepNumber)=iteration;
RunInfo.Forward.hResidual(CtrlVar.CurrentRunStepNumber)=r;


end


