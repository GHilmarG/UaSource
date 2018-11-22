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
MLC=BCs2MLC(MUA,BCs1);
Luv=MLC.ubvbL;
cuv=MLC.ubvbRhs;
Lh=MLC.hL;
ch=MLC.hRhs;

if numel(l1.ubvb)~=numel(cuv) ; l1.ubvb=zeros(numel(cuv),1) ; end
if numel(l1.h)~=numel(ch) ; l1.h=zeros(numel(ch),1) ; end
nlubvb=numel(l1.ubvb) ; 


[L,cuvh,luvh]=AssembleLuvh(Luv,Lh,cuv,ch,l1.ubvb,l1.h,MUA.Nnodes);
dl=luvh*0;


CtrlVar.uvhMatrixAssembly.ZeroFields=true; 
CtrlVar.uvhMatrixAssembly.Ronly=true;
[UserVar,RunInfo,R0,~]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1);
Fext0=R0;

iteration=0 ;
r=1e10; diffVector=zeros(CtrlVar.NRitmax,1); diffDu=1e10 ; diffDh=1e10; diffDlambda=1e10;  gamma=NaN;
RunInfo.Forward.Converged=0;
RunInfo.BackTrack.Converged=1 ; 

while true

    
    ResidualsCriteria=(~(r< CtrlVar.NLtol )) ...
        || iteration < CtrlVar.NRitmin;
    
    
    IncrementCriteria=(~(diffDu < CtrlVar.du && diffDh< CtrlVar.dh  )) ...
        || iteration < CtrlVar.NRitmin;
    
    
    
    if iteration > CtrlVar.NRitmax
        
        if CtrlVar.InfoLevelNonLinIt>=1
            fprintf(' SSTREAM(uvh) (time|dt)=(%g|%g): Maximum number of non-linear iterations reached. uvh iteration did not converge! \n',CtrlVar.time,CtrlVar.dt)
            fprintf(' Exiting uvh iteration after %-i iterations with r=%-g \n',iteration,r)
        end
        
        if CtrlVar.WriteRunInfoFile
            fprintf(RunInfo.File.fid,' SSTREAM(uvh) (time|dt)=(%g|%g): Maximum number of non-linear iterations reached. uvh iteration did not converge! \n',CtrlVar.time,CtrlVar.dt);
            fprintf(RunInfo.File.fid,' Exiting uvh iteration after %-i iterations with r=%-g \n',iteration,r);
        end
        
        RunInfo.Forward.Converged=0;
        break
    end
    
    if RunInfo.BackTrack.Converged==0
        if CtrlVar.InfoLevelNonLinIt>=1
            fprintf(' SSTREAM(uvh) (time|dt)=(%g|%g): Backtracting within non-linear iteration stagnated! \n Exiting non-lin iteraton with r=%-g, du=%-g and dh=%-g  after %-i iterations. \n',...
                CtrlVar.time,CtrlVar.dt,r,diffDu,diffDh,iteration) ;
        end
        
        if CtrlVar.WriteRunInfoFile
            fprintf(RunInfo.File.fid,' SSTREAM(uvh) (time|dt)=(%g|%g): Backtracting within non-linear iteration stagnated! \n Exiting non-lin iteraton with r=%-g, du=%-g and dh=%-g  after %-i iterations. \n',...
                CtrlVar.time,CtrlVar.dt,r,diffDu,diffDh,iteration) ;
        end
        
        RunInfo.Forward.Converged=0;
        break
    end
    
        
    switch lower(CtrlVar.uvhConvergenceCriteria)
        
        case 'residuals'
            
            if ~ResidualsCriteria
                
                tEnd=toc(tStart);
                if CtrlVar.InfoLevelNonLinIt>=1
                    fprintf(' SSTREAM(uvh) (time|dt)=(%g|%g): Converged to given residual tolerance of %-g with r=%-g in %-i iterations and in %-g  sec \n',...
                        CtrlVar.time,CtrlVar.dt,CtrlVar.NLtol,r,iteration,tEnd) ;
                end
                RunInfo.Forward.Converged=1;
                break
                
            end
            
        case 'increments'
            
            if ~IncrementCriteria
                
                tEnd=toc(tStart);
                if CtrlVar.InfoLevelNonLinIt>=1
                    fprintf(' SSTREAM(uvh) (time|dt)=(%g|%g): Converged to given increment tolerance of du=%g and dh=%g with r=%-g in %-i iterations and in %-g  sec \n',...
                        CtrlVar.time,CtrlVar.dt,CtrlVar.du,CtrlVar.dh,r,iteration,tEnd)
                end
                RunInfo.Forward.Converged=1;
                break
                
            end
            
        case 'residuals and increments'
            
            if ~ResidualsCriteria && ~IncrementCriteria
                
                tEnd=toc(tStart);
                if CtrlVar.InfoLevelNonLinIt>=1
                fprintf(CtrlVar.fidlog,' SSTREAM(uvh) (time|dt)=(%g|%g): Converged to given residual (r=%g) and increment tolerances (du=%g,dh=%g) with r=%-g, du=%-g and dh=%-g in %-i iterations and in %-g  sec \n',...
                    CtrlVar.time,CtrlVar.dt,CtrlVar.NLtol,CtrlVar.du,CtrlVar.dh,r,diffDu,diffDh,iteration,tEnd) ;
                end
                RunInfo.Forward.Converged=1;
                break
            end
            
         case 'residuals or increments'
            
            if ~ResidualsCriteria || ~IncrementCriteria
                
                tEnd=toc(tStart);       
                if CtrlVar.InfoLevelNonLinIt>=1
                  fprintf(CtrlVar.fidlog,' SSTREAM(uvh) (time|dt)=(%g|%g): Converged to given residual (r=%g) or increment tolerances (du=%g,dh=%g) with r=%-g, du=%-g and dh=%-g in %-i iterations and in %-g  sec \n',...
                    CtrlVar.time,CtrlVar.dt,CtrlVar.NLtol,CtrlVar.du,CtrlVar.dh,r,diffDu,diffDh,iteration,tEnd) ;
                end
                RunInfo.Forward.Converged=1;
                break
            end    
        otherwise
            fprintf(' CtrlVar.uvhConvergenceCriteria (%s) not set to a valid value.\n',CtrlVar.uvhConvergenceCriteria)
            error('parameter values incorrect')
    end



    iteration=iteration+1;
        
 
    %% Residuals , at gamma=0;
    gamma=0;
    [UserVar,RunInfo,r0,ruv0,rh0,rl0,R,K,frhs,grhs]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,Fext0);
    [duvh,dl]=solveKApe(K,L,frhs,grhs,[dub;dvb;dh],dl,CtrlVar);
    dub=duvh(1:MUA.Nnodes) ;  dvb=duvh(MUA.Nnodes+1:2*MUA.Nnodes); dh=duvh(2*MUA.Nnodes+1:end);
  
    %% calculate  residuals at full Newton step, i.e. at gamma=1
    gamma=1;
        
    [UserVar,RunInfo,r1,ruv1,rh1,rl1]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,Fext0);
    
    %% either accept full Newton step or do a line search

    
    
    [UserVar,RunInfo,gamma,r,ruv,rh,rl]=FindBestGamma2DuvhBacktrack(UserVar,RunInfo,CtrlVar,MUA,F0,F1,dub,dvb,dh,dl,L,luvh,cuvh,r0,r1,ruv1,rh1,rl1,Fext0);
    
    % [UserVar,RunInfo,r1Test,ruv1Test,rh1Test,rl1Test]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,Fext0);
    %
    % Is r1Test equal to r ?
    %
    
    iarm=RunInfo.BackTrack.iarm;
    infovector=RunInfo.BackTrack.Infovector;
    

    
    
    %% If desired, plot residual along search direction
    if CtrlVar.InfoLevelNonLinIt>=10 && CtrlVar.doplots==1
        nnn=12;
        gammaTestVector=zeros(nnn,1) ; rTestvector=zeros(nnn,1);
        Up=2.2;
        if gamma>0.7*Up ; Up=2*gamma; end
        parfor I=1:nnn
            gammaTest=Up*(I-1)/(nnn-1)+gamma/50;
            [~,~,rTest,~,~,~]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gammaTest,Fext0);
            gammaTestVector(I)=gammaTest ; rTestvector(I)=rTest;
        end
        
        gammaTestVector=[gammaTestVector(:);infovector(:,1)];
        rTestvector=[rTestvector(:);infovector(:,2)];
        
        [gammaTestVector,ind]=unique(gammaTestVector) ; rTestvector=rTestvector(ind) ;
        [gammaTestVector,ind]=sort(gammaTestVector) ; rTestvector=rTestvector(ind) ;
        
        slope=-2*r0;
        
        FigName='uvh iteration: line-search';
        fig=findobj(0,'name',FigName);
        if isempty(fig)
            fig=figure('name',FigName);
            fig.Position=[10,10,600,600] ;
        else
            fig=figure(fig);
            hold off
        end
        
        
        
        plot(gammaTestVector,rTestvector,'o-r') ; hold on ;
        plot(gamma,r,'Marker','h','MarkerEdgeColor','k','MarkerFaceColor','g')
        plot([gammaTestVector(1) gammaTestVector(2)],[rTestvector(1) rTestvector(1)+(gammaTestVector(2)-gammaTestVector(1))*slope],'g')
        
        title(sprintf('uvh iteration %-i,  iarm=%-i ',iteration,iarm)) ; xlabel(' \gamma ') ; ylabel('Residual')
        
        
        
        
        
        hold off
        %input('press return to continue')
    end
    
    
    %% calculate statistics on change in speed, thickness and Lagrange parameters
    D=mean(sqrt(F1.ub.*F1.ub+F1.vb.*F1.vb))+CtrlVar.SpeedZero;
    diffDu=full(max(abs(dub))+max(abs(dvb)))/D;        % sum of max change in du and dv normalized by mean speed
    diffDh=full(max(abs(dh))/mean(abs(F1.h)));            % max change in thickness divided by mean thickness
    diffDlambda=max(abs(dl))/mean(abs(luvh));
    diffVector(iteration)=r0;   % override last value, because it was just an (very accurate) estimate
    diffVector(iteration+1)=r;
    
    
    if isempty(diffDlambda)
        diffDlambda=0;
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
    [F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);
    CtrlVar.ResetThicknessToMinThickness=temp;
    
    if~isempty(Lh)
        BCsNormh=norm(ch-Lh*F1.h);
    else
        BCsNormh=0;
    end
    
    if~isempty(Luv)
        BCsNormuv=norm(cuv-Luv*[F1.ub;F1.vb]);
    else
        BCsNormuv=0;
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
        PlotForceResidualVectors('uvh',R,L,luvh,MUA.coordinates,CtrlVar) ; axis equal tight
    end
    
    if CtrlVar.InfoLevelNonLinIt>=1
        fprintf(...
            'NR-STREAM(uvh):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , ruv=%-14.7g , rh=%-14.7g , du=%-14.7g , dh=%-14.7g , dl=%-14.7g , BCsNormuv=%-g , BCsNormh=%-g  \n ',...
            iteration,iarm,gamma,r/r0,r0,r,ruv,rh,diffDu,diffDh,diffDlambda,BCsNormuv,BCsNormh);
        
    end
    
    
    
    if CtrlVar.WriteRunInfoFile
        
        fprintf(RunInfo.File.fid,...
            'NR-STREAM(uvh):%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , ruv=%-14.7g , rh=%-14.7g , du=%-14.7g , dh=%-14.7g , dl=%-14.7g , BCsNormuv=%-g , BCsNormh=%-g  \n ',...
            iteration,iarm,gamma,r/r0,r0,r,ruv,rh,diffDu,diffDh,diffDlambda,BCsNormuv,BCsNormh);
        
    end
    
    
    
    
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
    
    N=max([1,iteration-5]);
    
    [~,~,a1]=detrend_xt(log10(diffVector(N:iteration)),N:iteration);
    fprintf(CtrlVar.fidlog,' slope NR : %14.7g \n',a1);
    
    FigName='NR uvh implicit';
    fig=findobj(0,'name',FigName);
    if isempty(fig)
        fig=figure('name',FigName);
        fig.Position=[10,10,800,800] ;
    else
        fig=figure(fig);
        hold off
    end
    
    
    semilogy(0:iteration,diffVector(1:iteration+1),'x-r') ; title('NR uvh implicit') ; xlabel('Iteration') ; ylabel('Residual')
end

if ~isempty(L)
    BCerror=norm(L*[F1.ub;F1.vb;F1.h]-cuvh);
    if BCerror>0
        fprintf(CtrlVar.fidlog,'Norm of error satisfying Dirichlet BC=%14.7g  \n ',norm(L*[F1.ub;F1.vb;F1.h]-cuvh));
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

RunInfo.Forward.Iterations=iteration;
RunInfo.Forward.Residual=r;
RunInfo.Forward.IterationsTotal=RunInfo.Forward.IterationsTotal+RunInfo.Forward.Iterations; 

if CtrlVar.WriteRunInfoFile
    
    fprintf(RunInfo.File.fid,' --->  SSTREAM(uvh/%s) \t time=%15.5f \t dt=%-g \t r=%-g \t #it=% i \t CPUsec=%-g \n',...
        CtrlVar.uvhTimeSteppingMethod,CtrlVar.time,CtrlVar.dt,RunInfo.Forward.Residual,RunInfo.Forward.Iterations,tEnd) ;
    
end


end


