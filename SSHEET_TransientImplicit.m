function [UserVar,u1,v1,h1,s1,GF1,lambdah1,RunInfo]=SSHEET_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,BCs,dt,h1,h0,S,B,as0,ab0,as1,ab1,lambdah,AGlen,n,rho,rhow,g)

%  s and h are the initial estimates for s1 and h1
%  these are then updated, once convergent, I set s1=s and h1=h


nargoutchk(8,8)
narginchk(20,20)


[b0,s0,h0,~]=Calc_bs_From_hBS(CtrlVar,MUA,h0,S,B,rho,rhow);
[b1,s,h,~]=Calc_bs_From_hBS(CtrlVar,MUA,h1,S,B,rho,rhow);



a0=as0+ab0 ; a1=as1+ab1;
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


if any(h0<0) ; warning('MATLAB:uvh2DSSHEET',' thickness negative ') ; end

tStart=tic;

MLC=BCs2MLC(CtrlVar,MUA,BCs);
Lh=MLC.hL ; ch=MLC.hRhs ;

if numel(lambdah)~=numel(ch) ; lambdah=zeros(numel(ch),1) ; end

dlambdah=lambdah*0;
dh=h0*0;
iteration=0 ; Stagnated=0;
r=1e10; diffVector=zeros(CtrlVar.NRitmax,1) ; diffDh=1e10; diffDlambda=1e10;


while ((r> CtrlVar.NLtol || diffDh> CtrlVar.dh  || diffDlambda > CtrlVar.dl) &&  iteration <= CtrlVar.NRitmax && ~Stagnated)  || iteration < CtrlVar.NRitmin
    
    
    if (r> CtrlVar.NLtol) &&  (diffDh <  1e-6)
        fprintf(' (r> CtrlVar.NLtol) &&  (diffDh <  1e-6) ') 
    end
    
    iteration=iteration+1;
    
    %fprintf('It:%-i numel(dlambdah)=%-i \t numel(lambdah)=%-i \n ',iteration,numel(dlambdah),numel(lambdah))
    
    [R,K,F,T]=MatrixAssemblySSHEETtransient2HD(CtrlVar,MUA,AGlen,n,rho,g,s0,b0,s,b1,a0,a1,dt);
    
    if iteration==1
        F0=F ;  % There is a potential issue here which is that F0 is zero if the accumulation
                % and grad q is everywhere zero. However, if this happens dh/dt will also automatically be
                % zero. 
                % F0 is used as a normalisation factor when calculating the residual,
                % do not change this normalisation factor in the course of the iteration
    end

    gamma=0;
    
    if ~isempty(Lh)
        r0=ResidualCostFunction(R+Lh'*(lambdah+gamma*dlambdah),[],F0,MUA.Nnodes);
    else
        r0=ResidualCostFunction(R,[],F0,MUA.Nnodes);
    end
    
    %% solve the linear system
    if ~isempty(Lh)
        frhs=-R-Lh'*lambdah;
        grhs=ch-Lh*h;
    else
        frhs=-R;
        grhs=[];
    end
    
    
    [dh,dlambdah]=solveKApe(K,Lh,frhs,grhs,dh,dlambdah,CtrlVar);
    
    
    if any(isnan(dh)) ; save TestSave  ;
        fprintf(CtrlVar.fidlog,'error: NaN in solution of implicit system \n') ;
        error(' NaN in solution of implicit h system ' ) ;
    end
    
    
    %% calculate  residuals at full Newton step
    
    func=@(gamma) CalcCostFunctionSSHEET(CtrlVar,gamma,dh,MUA,AGlen,n,rho,g,s0,b0,s,b1,a0,a1,dt,Lh,lambdah,dlambdah,F0,ch);
               
    gamma=1;
    r1=func(gamma);
    
    
    
    %% either accept full Newton step or do a line search
    %[r,gamma,infovector,iarm,BacktrackInfo]=FindBestGammaBacktrackSSHEET(CtrlVar,r0,r1,F0,MUA,AGlen,n,rho,g,s0,b0,s,b1,a0,a1,dh,dt,Lh,lambdah,dlambdah);
    
    %fprintf(' r=%-g \t gamma=%-g \n ',r,gamma)
    %F=@(q,u,v) func(C0-q*dJdC,AGlen,u,v); nOut=9; listInF=[1 2] ; listOutF=[7 8];
    Slope0=-2*r0 ;  % using the inner product def
    [gamma,r,BackTrackInfo]=BackTracking(Slope0,1,r0,r1,func,CtrlVar);
    
   
    
    infovector=BackTrackInfo.InfoVector;
    iarm=BackTrackInfo.nBackTrackSteps;
    
    %fprintf(' r=%-g \t gamma=%-g \n ',r,gamma)
    
    if BackTrackInfo.converged==0
        Stagnated=1;
    end
    
    
    
    
    %% If desired, plot residual along search direction
    if CtrlVar.InfoLevelNonLinIt>=10 && CtrlVar.doplots==1
        nnn=12;
        gammaTestVector=zeros(nnn,1) ; rTestvector=zeros(nnn,1);
        Up=2.2;
        if gamma>0.7*Up ; Up=2*gamma; end
        parfor I=1:nnn
            gammaTest=Up*(I-1)/(nnn-1)+gamma/50;
            rTest=CalcCostFunctionSSHEET(CtrlVar,gammaTest,dh,MUA,AGlen,n,rho,g,s0,b0,s,b1,a0,a1,dt,Lh,lambdah,dlambdah,F0,ch);
            gammaTestVector(I)=gammaTest ; rTestvector(I)=rTest;
        end
        gammaTestVector=[gammaTestVector(:);infovector(:,1)];
        rTestvector=[rTestvector(:);infovector(:,2)];
        [gammaTestVector,ind]=unique(gammaTestVector) ; rTestvector=rTestvector(ind) ;
        [gammaTestVector,ind]=sort(gammaTestVector) ; rTestvector=rTestvector(ind) ;
        
        
        slope=-2*r0;
        figure ; plot(gammaTestVector,rTestvector,'o-r') ; hold on ;
        plot(gamma,r,'Marker','h','MarkerEdgeColor','k','MarkerFaceColor','g')
        plot([gammaTestVector(1) gammaTestVector(2)],[rTestvector(1) rTestvector(1)+(gammaTestVector(2)-gammaTestVector(1))*slope],'g')
        
        title(sprintf('SSHEET h iteration %-i,  iarm=%-i ',iteration,iarm)) ; xlabel(' \gamma ') ; ylabel('Residual')
        
        hold off
        %input('press return to continue')
    end
    
    
   
    
    %% update variables
    h=h+gamma*dh ; lambdah=lambdah+gamma*dlambdah;
    
    temp=CtrlVar.ResetThicknessToMinThickness;
    if ~CtrlVar.ResetThicknessInNonLinLoop
        CtrlVar.ResetThicknessToMinThickness=0;
    end
    [b1,s,h,~]=Calc_bs_From_hBS(CtrlVar,MUA,h,S,B,rho,rhow);
    CtrlVar.ResetThicknessToMinThickness=temp;
    %[b1,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
    
    
     %% calculate statistics on change in speed, thickness and Lagrange parameters
    
    isThickPos=h>=CtrlVar.ThickMin; 
    diffDh=norm(dh(isThickPos))/norm(h(isThickPos))  ;       % norm in change of the solution during this iteration, where thick positive.
    
    % diffDh=gamma*full(max(abs(dh))/max(abs(h0)));            % max change in thickness divided by mean thickness
    
    diffDlambda=gamma*full(max(abs(dlambdah))/max(abs(lambdah)));
    if~isempty(Lh)
        BCsNorm=norm(ch-Lh*h);
    else
        BCsNorm=0;
    end
    %fprintf(' BCsNorm=%-g \n ',BCsNorm)
    diffVector(iteration)=r0;   % override last value, because it was just an (very accurate) estimate
    diffVector(iteration+1)=r;
    
    if isempty(diffDlambda)
        diffDlambda=0;
    end
    
    
    %% plot and print info
    if CtrlVar.InfoLevelNonLinIt>=100  && CtrlVar.doplots==1
        PlotForceResidualVectors('h-only',R,Lh,lambdah,MUA.coordinates,CtrlVar) ; axis equal tight
    end
    
    
    if CtrlVar.InfoLevelNonLinIt>=1
        fprintf(CtrlVar.fidlog,'NR-SSHEET(h):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , dh=%-14.7g , dl=%-14.7g , BCsNorm=%-14.7g\n ',...
            iteration,iarm,gamma,r/r0,r0,r,diffDh,diffDlambda,BCsNorm);
    end
    
    
end

%% return calculated values at the end of the time step
h1=h ; lambdah1=lambdah;
[~,s1,h1,GF1]=Calc_bs_From_hBS(CtrlVar,MUA,h1,S,B,rho,rhow);
%[b1,s1,h1]=Calc_bs_From_hBS(h1,S,B,rho,rhow,CtrlVar,MUA.coordinates);
[u1,v1]=uvSSHEET(CtrlVar,MUA,BCs,AGlen,n,rho,g,s1,h1);

tEnd=toc(tStart);

%% print/plot some info
if CtrlVar.InfoLevelNonLinIt>=1 && r < CtrlVar.NLtol 
    fprintf(CtrlVar.fidlog,' SSHEET(h) converged to given tolerance of %-g with r=%-g in %-i iterations and in %-g  sec \n',CtrlVar.NLtol,r,iteration,tEnd) ;
end

if CtrlVar.InfoLevelNonLinIt>=10 && iteration >= 2 && CtrlVar.doplots==1
    
    N=max([1,iteration-5]);
    
    [~,~,a1]=detrend_xt(log10(diffVector(N:iteration)),N:iteration);
    fprintf(CtrlVar.fidlog,' slope NR : %14.7g \n',a1);
    figure; semilogy(0:iteration,diffVector(1:iteration+1),'x-r') ; title('NR SSHEET h implicit') ; xlabel('Iteration') ; ylabel('Residual')
end

%     if ~isempty(Lh) &&   CtrlVar.InfoLevelNonLinIt>=0
%         fprintf(CtrlVar.fidlog,' final error in satisfying Dirichlet BC %14.7g  \n ',norm(Lh*h-Lhrhs));
%     end

RunInfo.Forward.Converged=1;
if r>CtrlVar.NLtol
    fprintf(CtrlVar.fidlog,' SSHEET  did not converge to given tolerance of %-g with r=%-g in %-i iterations and in %-g  sec \n',CtrlVar.NLtol,r,iteration,tEnd);
    RunInfo.Forward.Converged=0;
end

if iteration > CtrlVar.NRitmax
    fprintf(CtrlVar.fidlog,'Warning: maximum number of NRh iterations %-i reached \n',CtrlVar.NRitmax);
    warning('uvh2DSSHEET:MaxIterationReached','uvh2DSSHEET exits because maximum number of iterations %-i reached \n',CtrlVar.NRitmax)
end

RunInfo.Forward.Iterations=iteration;   
RunInfo.Forward.IterationsTotal=RunInfo.Forward.IterationsTotal+RunInfo.Forward.Iterations; 
RunInfo.Forward.Residual=r;


end


