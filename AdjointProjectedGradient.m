function   [UserVar,Cest,AGlenEst,Info,ub,vb,ud,vd,xAdjoint,yAdjoint,gammaAdjoint]=AdjointProjectedGradient(...
    UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info)

%
%                [Cest,AGlenEst,Info,ub,vb,ud,vd,xAdjoint,yAdjoint,gammaAdjoint]=AdjointProjectedGradient(...
%                 UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info);

% Minimisation using the projected (conjugated) gradient method.

if CtrlVar.InfoLevelAdjoint>=10
    fprintf(' minimisation method AdjointProjectedGradient \n ')
end


AGlenEst=InvStartValues.AGlen;
Cest=InvStartValues.C;
n=InvStartValues.n;
m=InvStartValues.m;

upA=zeros(length(AGlenEst),1)+CtrlVar.AGlenmax;  lowA=zeros(length(AGlenEst),1)+CtrlVar.AGlenmin;
upC=zeros(length(Cest),1)+CtrlVar.Cmax;  lowC=zeros(length(Cest),1)+CtrlVar.Cmin;

if norm(Cest - kk_proj(Cest,upC,lowC)) > 0
    disp(' initial C iterate not feasible ');
    Cest=kk_proj(Cest,upC,lowC);
end

if norm(AGlenEst - kk_proj(AGlenEst,upA,lowA)) > 0
    disp(' initial AGlen iterate not feasible ');
    AGlenEst=kk_proj(AGlenEst,upA,lowA);
end


dJdAdescent=AGlenEst*0;
dJdCdescent=Cest*0;

nIt=CtrlVar.MaxAdjointIterations;

if isempty(Info.JoptVector)
    Info.JoptVector=zeros(nIt,7)+NaN; iJ=0;
else
    iJ=size(Info.JoptVector,1);
    Info.JoptVector=[Info.JoptVector;zeros(nIt,7)+NaN];
end

gammaAdjoint=0;

for iteration=1:nIt
    
    % At the beginning of the iteration the best estimates for AGlen and C are AGlenEst and Cest
    % These values are the starting values of the iteration, ie C0 and AGlen0
    % In the line search I modify C0 and AGlen0 and these modified values are Ctest and AGlentest
    % Once the line search is finised, AGlenEst and Cest are updated
    
    %if CtrlVar.doplots==1 && mod(iteration,6)==0 ; close all ; end
    
    iA=strfind(CtrlVar.AdjointGrad,'A'); iC=strfind(CtrlVar.AdjointGrad,'C'); isAgrad=~isempty(iA); isCgrad=~isempty(iC);
    
    
    if isAgrad && isCgrad
        CtrlVar.AdjointConjugatedGradients=0 ;
        if mod(iteration,2)
            isAgrad=1 ; isCgrad=0 ;
        else
            isAgrad=0 ; isCgrad=1  ;
        end
    end
    
    
    if isAgrad ; fprintf(' A gradient optimisation step \n ') ; end
    if isCgrad ; fprintf(' C gradient optimisation step \n ') ; end
    
    
    ig=0; fVector=zeros(100,6)+NaN; gamma_Vector=zeros(100,1)+NaN;
    C0=Cest;
    AGlen0=AGlenEst;
    
    [J0,Idata0,IRegC0,IRegAGlen0,IBarrierC0,IBarrierAGlen0,ub,vb,ud,vd,l,dIdu0,Kuv,Ruv,RunInfo]=...
        CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen0,C0,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
    
    
    
    JMin=J0 ; IdataMin=Idata0 ;
    IRegCmin=IRegC0; IRegAGlenmin=IRegAGlen0;
    IBarrierCmin=IBarrierC0 ; IBarrierAGlenmin=IBarrierAGlen0;
    
    if iJ==0
        iJ=iJ+1;  
        Info.JoptVector(iJ,1)=J0; 
        Info.JoptVector(iJ,2)=Idata0;
        Info.JoptVector(iJ,3)=IRegC0; 
        Info.JoptVector(iJ,4)=IRegAGlen0;
        Info.JoptVector(iJ,5)=IBarrierC0; 
        Info.JoptVector(iJ,6)=IBarrierAGlen0;
        Info.JoptVector(iJ,7)=NaN;
        Info.AdjointGrad{iJ}=CtrlVar.AdjointGrad;
    end
    
    ig=ig+1; fVector(ig,1)=J0; fVector(ig,2)=Idata0 ;
    fVector(ig,3)=IRegC0;  fVector(ig,4)=IRegAGlen0;
    fVector(ig,5)=IBarrierC0 ; fVector(ig,6)=IBarrierAGlen0 ;
    gamma_Vector(ig)= 0;
    
    %% Gradient calculated using brute force method (for verification purposes)
    if CtrlVar.CalcBruteForceGradient
        
        [dJdCBruteForce,dJdAGlenBruteForce]=...
            CalcBruteForceGradient(CtrlVar,MUA,s,S,B,h,ub,vb,AGlen0,C0,Luv,Luvrhs,ubvbLambda,n,m,alpha,rho,rhow,g,...
            uMeas,vMeas,wMeas,C_prior,AGlen_prior,Cd,CAGlen,CC,GF);
        
        if  isCgrad &&  CtrlVar.InfoLevelAdjoint>=10 && CtrlVar.doplots==1
            figure
            PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,dJdCBruteForce); colorbar ; axis equal tight
            title('Brute force dJdC') ; xlabel('x') ; ylabel('y')
        end
        
        if  isAgrad &&  CtrlVar.InfoLevelAdjoint>=10 && CtrlVar.doplots==1
            PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates,dJdAGlenBruteForce);
            title('Brute force dJdAGlen') ; xlabel('x') ; ylabel('y')
        end
        
    end
    
    %% Adjoint method
    
    [UserVar,dJdC,dJdAGlen,ub,vb,ud,vd,xAdjoint,yAdjoint,dIdCreg,dIdAGlenreg,dIdCdata,dIdAGlendata,dIdCbarrier,dIdAGlenbarrier,lambdaAdjoint]=...
        AdjointGradientNR2d(...
        UserVar,CtrlVar,MUA,BCs,BCsAdjoint,s,b,h,S,B,ub,vb,ud,vd,l,AGlen0,C0,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
        
    indA=AGlen0 == upA | AGlen0 ==lowA ;
    indC=C0 == upC | C0 ==lowC ;
    fprintf(' fraction of active AGlen constrains %-g \n ',sum(indA)/length(AGlen0))
    fprintf('     fraction of active C constrains %-g \n ',sum(indC)/length(C0))
    
    %%
    % Project the gradient: if at the boundary and corresponding elements of (negative) gradient pointing ouside,
    % set to zero
    % (only important in connection with conjugated gradient option)
    
    IlowerC=(C0<=lowC) & (dJdC>0);
    IupperC=(C0>=upC)  & (dJdC<0);
    
    fprintf(' fraction of lower and upper eps active AGlen constrains %-g \t %-g \n ',sum(IlowerC)/length(C0),sum(IupperC)/length(C0))
    
    dJdC(IlowerC)=0; % at lower limit and (negative) gradient points out of feasible domain
    dJdC(IupperC)=0;
    
    IlowerA=(AGlen0<=lowA) & (dJdAGlen>0) ;
    IupperA=(AGlen0>=upA)  & (dJdAGlen<0) ;
    
    fprintf(' fraction of lower and upper eps active AGlen constrains %-g \t %-g \n ',sum(IlowerA)/length(AGlen0),sum(IupperA)/length(AGlen0))
    
    dJdAGlen(IlowerA)=0;
    dJdAGlen(IupperA)=0;
    
    if ~isCgrad ;  dJdCsearch=zeros(length(Cest),1)     ; dJdC=zeros(length(Cest),1)     ; end
    if ~isAgrad ;  dJdAGlensearch=zeros(length(AGlenEst),1) ; dJdAGlen=zeros(length(AGlenEst),1) ; end
    
    if CtrlVar.CalcBruteForceGradient
        save AdjointBruteForceGradients  MUA dJdAGlen dJdC dJdCBruteForce dJdAGlenBruteForce dIdAGlenreg dIdCreg CtrlVar
        %           error(' no error just feel like stopping ' )
    end
    
    %       save AdjointGradients coordinates connectivity dJdAGlen dJdC dIdAGlenreg dIdCreg CtrlVar
    %       error(' no error just feel like stopping ' )
    
    %%
    
    
    %%  Testing if the direction dIdatadC is a desent direction, and correcting the magnitude of the gradient
    % I=I(C0+gamma dIdCest)  => |dIdgamma|=dIdC'*dIdCest
    % define kappa as dIdC=kappa dIdCest (assuming it is in same direction)
    % if I calculate dIdgamma numerically I find: dIdgamma=dIdC'*dIdCest=kappa*dIdCest'*dIdCest and that
    % kappa=dIdgamma/(dIdCest'*dIdCest)
    % When I calculate dIdgamma numerically I need a small perturbation, if the size of that step is gamma_eps then
    % the step size after scaling is gamma_eps/kappa.
    % After having found the min I can test if that step size is much smaller than the gamma that gives the minumum value
    % If not then I used too large step size in calculating the slope.
    %
    
    %dIdC=zeros(length(C),1) ; dIdA=zeros(length(AGlen),1) ;
    
    if isCgrad
        %% rescale C gradient
        gamma_eps=1e-4*norm(C0)/norm(dJdC); % small perturbation to C0
        if CtrlVar.RescaleAdjointGradient
            
            gamma=gamma_eps;
            Ctest=kk_proj(C0-gamma*dJdC,upC,lowC);  % changing C in steepest decent direction
            
            [Jeps,Idataeps,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,l,dIdu,Kuv,Ruv,RunInfo]=...
                CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen0,Ctest,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
            
            %[ub,vb,ubvbLambda]=SSTREAM2dNR(CtrlVar,MUA,s,S,B,h,ub,vb,AGlen0,Ctest,Luv,Luvrhs,ubvbLambda,n,m,alpha,rho,rhow,g);
            %[Jeps,Idataeps]=MisfitFunction(CtrlVar,MUA,us,vs,ws,sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,Cest,C_prior,AGlenEst,AGlen_prior,Cd,CAGlen,CC,GF);
            
            
            dJdgamma=(Jeps-J0)/gamma_eps;
            %dIdatadgamma=(Idataeps-Idata0)/gamma_eps;
            %             if dIdatadgamma >=0
            %                 warning('MinimisationSimple:desent','not a desent direction, i.e. going in the direction dIdatadC does not decrease the Idata cost function. ')
            %                 save TestSave
            %                 error('afsd')
            %             end
            
            scale=abs(dJdgamma)/(dJdC'*dJdC); % this involves an approximation assuming that dIdCest is in same direction as dIdC
            if CtrlVar.InfoLevelAdjoint>=10
                fprintf(' gradient scaling factor is %-g \n ',scale)
            end
            dJdC=scale*dJdC;
            
        end
        
        dJdCdescentlast=dJdCdescent;
        %dIdCdescent=-dIdatadC+CtrlVar.isRegC*dIdCreg;
        dJdCdescent=dJdC;
        
        %
        if CtrlVar.AdjointConjugatedGradients==1 && iteration>1
            dJdCsearchLast=dJdCsearch;
            [dJdCsearch,ConjGradAngle,teta]=NewConjugatedGrad(dJdCdescent,dJdCdescentlast,dJdCsearchLast,CtrlVar);
            
        else
            dJdCsearch=-dJdC;
        end
        
        dJdCsearch(AGlen0<=lowA & dJdCsearch<0)=0;  % project element of gradients of active set when point out of region
        dJdCsearch(AGlen0>=upA  & dJdCsearch>0)=0;
        
        
        
    end
    
    
    if isAgrad
        gamma_eps=1e-4*norm(AGlen0)/norm(dJdAGlen); % small perturbation to AGlen0
        if CtrlVar.RescaleAdjointGradient
            %% rescale AGlen gradient
            
            gamma=gamma_eps;
            AGlentest=kk_proj(AGlen0-gamma*dJdAGlen,upA,lowA);  % changing AGlen in the steepest decent direction
            
            [Jeps,Idataeps,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,l,dIdu,Kuv,Ruv,RunInfo]=...
                CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlentest,C0,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
            
            % [ub,vb,ubvbLambda]=SSTREAM2dNR(CtrlVar,MUA,s,S,B,h,ub,vb,AGlentest,C0,Luv,Luvrhs,ubvbLambda,n,m,alpha,rho,rhow,g);
            %[Jeps,Idataeps]=MisfitFunction(CtrlVar,MUA,us,vs,ws,sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,C0,C_prior,AGlentest,AGlen_prior,Cd,CAGlen,CC,GF);
            
            
            dJdgamma=(Jeps-J0)/gamma_eps;
            if dJdgamma >=0
                warning('MinimisationSimple:desent','not a desent direction, i.e. going in the direction dIdatadA does not decrease the Idata cost function. ')
            end
            
            scale=abs(dJdgamma)/(dJdAGlen'*dJdAGlen); % this involves an approximation
            fprintf(' scale for A is %-g \n ',scale)
            dJdAGlen=scale*dJdAGlen;
            
            
        end
        
        
        dJdAdescentlast=dJdAdescent;
        dJdAdescent=dJdAGlen;
        
        if CtrlVar.AdjointConjugatedGradients==1 && iteration>1
            dJdAsearchLast=dJdAGlensearch;
            [dJdAGlensearch,ConjGradAngle,teta]=NewConjugatedGrad(dJdAdescent,dJdAdescentlast,dJdAsearchLast,CtrlVar);
            
        else
            dJdAGlensearch=-dJdAdescent; %  dIdatadA-CtrlVar.isRegAGlen*dIdAreg;
        end
        dJdAGlensearch(AGlen0<=lowA & dJdAGlensearch <0)=0;  % project element of gradients of active set when point out of region
        dJdAGlensearch(AGlen0>=upA & dJdAGlensearch > 0)=0;
        
    end
    
    % Now dJdCsearch and dJdAGlensearch have bee determined
    
    
    
    %% simple line search
    
    if strcmpi(CtrlVar.AdjointGrad,'A')
        CtrlVar.AdjointGrad='A';
        Slope0=-dJdAGlen'*dJdAGlensearch;
    elseif strcmpi(CtrlVar.AdjointGrad,'C')
        Slope0=-dJdC'*dJdCsearch;
    else
        fprintf('CtrlVar.AdjointGrad=%s\',CtrlVar.AdjointGrad)
        error('Ua:AdjointProjectedGradient',' CtrlVar.AdjointGrad must be either equal to A or C')
    end
    
    
    fprintf('Slope0=%-g \n',Slope0)
    
    
    if iteration==1
        if ~isempty(CtrlVar.AdjointInitialSearchStepSize) && ~isnan(CtrlVar.AdjointInitialSearchStepSize) && isfinite(CtrlVar.AdjointInitialSearchStepSize)
            gamma_MinEstimate=CtrlVar.AdjointInitialSearchStepSize;
        elseif ~isempty(InvStartValues.InitialSearchStepSize) && ~isnan(InvStartValues.InitialSearchStepSize) && isfinite(InvStartValues.InitialSearchStepSize)
            gamma_MinEstimate=InvStartValues.InitialSearchStepSize ;
        else
            kappa=0.1 ; 
            %gamma_MinEstimate=kappa*J0/Slope0;
            gamma_MinEstimate=1;
        end
    else
        gamma_MinEstimate=gammaAdjoint; % last gammaAdjoint used as initial guess for step sized
    end
    
    
    gammaAdjoint=0 ; % now set gammaAdjoint to zero since at the start of line search this is where the current min is found (ie at gamma=0)
    
    if gamma_MinEstimate==0 ; gamma_MinEstimate=1 ; end
    % close to zero f=f0-Slope0*gamma,  aiming at reduction by kappa: kappa f0=f0-Slope gamma => gamma= f0 (1-kappa)/Slope
    %if iteration==1
    fprintf('\n ++++Inv. it. (start): %-i \t J/J0=%-g  \t J=%-g  \t J0=%-g \t Idata=%-g \t IRegC=%-g \t IBarrierC=%-g \t IRegAGlen=%-g \t IBarrierAGlen=%-g  \n \n ',...
        iteration+Info.InverseIterations,JMin/J0,JMin,J0,IdataMin,IRegCmin,IBarrierCmin,IRegAGlenmin,IBarrierAGlenmin)
    
    %
    %    fprintf('\n +++++++++++ At start of line search:  \t J0=%-g \t Idata=%-g \t RegC=%-g \t BarrierC=%-g \t RegAGlen=%-g \t BarrierAGlen=%-g \t gamma=%-g \n \n',...
    %        J0,Idata0,IRegC0,IBarrierC0,IRegAGlen0,IBarrierAGlen0,gamma_MinEstimate)
    %end
    
    gamma_a=0 ; Ja=J0 ; gamma_c=gamma_MinEstimate ; gamma_Eps=gamma_MinEstimate/1000;
    
    RunInfo.converged=0; iNR=0;
    while RunInfo.converged==0  && iNR<=5
        % possibly the change in C/AGlen is too large for the non-linear solver to converge
        % so I allow for a drastic reduciton in step size if needed
        iNR=iNR+1;
        gamma=gamma_c/10^(2*(iNR-1));
        gamma_b=gamma_c/2;
        
        Ctest=kk_proj(C0+gamma*dJdCsearch,upC,lowC);
        AGlentest=kk_proj(AGlen0+gamma*dJdAGlensearch,upA,lowA);
        
        if CtrlVar.InfoLevelAdjoint>10 
            fprintf('Solving forward problem using a line-search stepsize g=%g \n',gamma)
        end
        
            
        [UserVar,ub,vb,ud,vd,l,Kuv,Ruv,RunInfo,ubvbL]=...
            uv(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlentest,Ctest,n,m,alpha,rho,rhow,g,GF);
        
        if CtrlVar.InfoLevelAdjoint>10  && ~RunInfo.converged
           fprintf('Forward step did not converge. Will reduce line-search step size.\n') 
        end
        
    end
    
    if RunInfo.converged==0
        error('SSTREAM2d did not converge')
    end
    
    
    [Jc,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,l,dIdu,Kuv,Ruv,RunInfo]=...
        CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlentest,Ctest,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
    
    Jvalue=Jc;
    ig=ig+1; fVector(ig,1)=Jc; fVector(ig,2)=Idata ; fVector(ig,3)=IRegC;
    fVector(ig,4)=IRegAGlen; gamma_Vector(ig)= gamma;
    fVector(ig,5)=IBarrierC ; fVector(ig,6)=IBarrierAGlen ;
    iarm=0 ; iarmmax=CtrlVar.AdjointMaxLineSearchIterations;
    TargetRatio=0.5;
    target=TargetRatio*J0;
    
    
    
    if  Jc<=target
        JMin=Jc ; IdataMin=Idata ; IRegCmin=IRegC;  IRegAGlenmin=IRegAGlen;
        IBarrierCmin=IBarrierC ; IBarrierAGlenmin=IBarrierAGlen;
        gammaAdjoint=gamma;
        if CtrlVar.InfoLevelAdjoint>=10
            fprintf(' Initial step accepted. Misfit reduced from %-g to  %-g \t ratio fMin/f0=%-g \n ',J0,JMin,JMin/J0)
        end
    else
        if CtrlVar.InfoLevelAdjoint>=10
            fprintf(' Initial step not accepted entering line search. f0=%-g \t fc=%-g \t ratio fc/f0=%-g \n ',J0,Jc,Jc/J0)
        end
        
        if Jc<J0
            JMin=Jc ; IdataMin=Idata ; IRegCmin=IRegC;  IRegAGlenmin=IRegAGlen;
            IBarrierCmin=IBarrierC ; IBarrierAGlenmin=IBarrierAGlen;
            gammaAdjoint=gamma;
        end
        
        iFminTry=1; iFminTryMax=3;
        % I continue linesearch if any of these conditions are fullfilled
        
        while JMin > target && iarm <= iarmmax && iFminTry <= iFminTryMax
            iarm=iarm+1;
            iFminTry=iFminTry+1;
            
            gamma=gamma_b;
            Ctest=kk_proj(C0+gamma*dJdCsearch,upC,lowC);
            AGlentest=kk_proj(AGlen0+gamma*dJdAGlensearch,upA,lowA);
            
            Jvaluelast=Jvalue;
            
            % [ub,vb]=SSTREAM2dNR(CtrlVar,MUA,s,S,B,h,ub,vb,AGlentest,Ctest,Luv,Luvrhs,ubvbLambda,n,m,alpha,rho,rhow,g);
            % [Jb,Idata,IRegC,IRegAGlen,dIduv,IBarrierC,IBarrierAGlen]=MisfitFunction(UserVar,CtrlVar,MUA,ub,vb,ud,vd,AGlen,C,Priors,Meas);
            
            [Jb,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,l,dIdu,Kuv,Ruv,RunInfo]=...
                CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlentest,Ctest,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
            
            
            ig=ig+1; fVector(ig,1)=Jb; fVector(ig,2)=Idata ;
            fVector(ig,3)=IRegC;  fVector(ig,4)=IRegAGlen;
            fVector(ig,5)=IBarrierC ; fVector(ig,6)=IBarrierAGlen ;
            gamma_Vector(ig)= gamma; Jvalue=Jb;
            
            if Jvalue < Jvaluelast || JMin >= J0 
                iFminTry=0;
            end   % reset iFminTry if function value is decreasing or if fMin is still larger than f0
            
            if Jb < JMin 
                JMin=Jb ;  IdataMin=Idata ; IRegCmin=IRegC; IRegAGlenmin=IRegAGlen;
                IBarrierCmin=IBarrierC ; IBarrierAGlenmin=IBarrierAGlen;
                gammaAdjoint=gamma_b ;
            end
            
            if CtrlVar.InfoLevelAdjoint>=10
                fprintf('Line-search step # %-i. f(a)=%-10.5g \t f(b)=%-10.5g \t f(c)=%-10.5g \t fmin=%-10.5g  \t fmin/ft=%-10.5g \t fmin/f0=%-g \n ',...
                    iarm,Ja,Jb,Jc,JMin,JMin/target,JMin/J0)
                fprintf('                         a=%-10.5g  \t   b=%-10.5g  \t    c=%-10.5g   \t    g=%-10.5g \n ',gamma_a,gamma_b,gamma_c,gammaAdjoint)
            end
            
            
            [gamma_Test,ParStatus ] = parabolamin(gamma_a,gamma_b,gamma_c,Ja,Jb,Jc);
            
            if ParStatus==1  % parabolic fit did not give an acceptable minimum
                if Jc < Jb && Jb < Ja      % decreasing values
                    gamma_Test=10*gamma_c;    % extrapolate
                elseif Jc > Jb && Jb > Ja  % increasing values
                    gamma_Test=(gamma_a+gamma_b)/10;
                else                       %
                    gamma_Test=(gamma_a+gamma_b)/3;
                end
                
                if isinf(Ja+Jb+Jc)
                    fprintf(' function values are Inf. Try to increase step size \n')
                    gamma_Test=10*gamma_c;    % extrapolate
                end
                
            end
            
            [gamma_Dist,imin]=min(abs(gamma_Test-gamma_Vector(2:ig-1)));
            if gamma_Dist < gamma_Eps
                if JMin>0.999*J0
                    gamma_Eps=gamma_Eps/100;
                else
                    if CtrlVar.InfoLevelAdjoint>=10
                        fprintf(' breaking out of linesearch because new value (%-g) so close to a prevous one (%-g) \n',gamma_Test,gamma_Vector(1+imin))
                    end
                    iarm=iarm-1;
                    break
                end
            end
            
            if gamma_Test<0 ; gamma_Test=(gamma_a+gamma_b)/10; end
            %if gamma_Min > 0.8*gamma_b ; gamma_Min=0.8*gamma_b ; elseif gamma_Min < 0.1*gamma_b ; gamma_Min=0.1*gamma_b; end
            
            if gamma_Test < gamma_b
                gamma_c=gamma_b ; Jc=Jb; gamma_b=gamma_Test ; Ja=J0 ; gamma_a=0;
            elseif gamma_Test >= gamma_b && gamma_Test< gamma_c
                gamma_a=gamma_b ; Ja=Jb ; gamma_b=gamma_Test ;
            else  % extrapolation step
                fprintf(' extrapolation step ')
                if gamma_Test> 5*gamma_c ; gamma_Test=5*gamma_c ; end
                gamma_a=gamma_b ; Ja=Jb ; gamma_b=gamma_Test ; gamma_c=5*gamma_Test;
                
                gamma=gamma_c;
                Ctest=kk_proj(C0+gamma*dJdCsearch,upC,lowC);
                AGlentest=kk_proj(AGlen0+gamma*dJdAGlensearch,upA,lowA);
                
                %[ub,vb]=SSTREAM2dNR(CtrlVar,MUA,s,S,B,h,ub,vb,AGlentest,Ctest,Luv,Luvrhs,ubvbLambda,n,m,alpha,rho,rhow,g);
                %[Jc,Idata,IRegC,IRegAGlen,dIduv,IBarrierC,IBarrierAGlen]=MisfitFunction(UserVar,CtrlVar,MUA,ub,vb,ud,vd,AGlentest,Ctest,Priors,Meas);
                
                [Jc,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,l,dIdu,Kuv,Ruv,RunInfo]=...
                    CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlentest,Ctest,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
                
                Jvaluelast=Jvalue ; Jvalue=Jc; if Jvalue < Jvaluelast  ; iFminTry=0 ; end
                if Jc <JMin ;
                    JMin=Jc ;
                    IdataMin=Idata ; IRegCmin=IRegC;  IRegAGlenmin=IRegAGlen;
                    IBarrierCmin=IBarrierC ; IBarrierAGlenmin=IBarrierAGlen;
                    gammaAdjoint=gamma_c ;
                end
                ig=ig+1; fVector(ig,1)=Jc; fVector(ig,2)=Idata ;
                fVector(ig,3)=IRegC;  fVector(ig,4)=IRegAGlen;
                fVector(ig,5)=IBarrierC ; fVector(ig,6)=IBarrierAGlen ;
                gamma_Vector(ig)= gamma;
                
            end
        end
        if ~(iarm<=iarmmax) ; fprintf(' Line search exit due to maximum number of iterations reached. \n') ; end
        if ~(iFminTry<=iFminTryMax) ; fprintf(' Exiting line search because no improvement found for the last %-i iterations. \n',iFminTry) ; end
        if ~(JMin >  target); fprintf(' Set target ratio-reduction of %-g reached in line search. \n',TargetRatio) ; end
        
    end
    
%% plotting
    if CtrlVar.InfoLevelAdjoint>=100 && CtrlVar.doplots==1;
        
        nLS=5;
        
        gamma_StepVector=gammaAdjoint*[-1/3,1/3,2/3,4/3,5/3];
        Jplot=zeros(nLS,1) ; Idata=zeros(nLS,1) ; IRegC=zeros(nLS,1); IRegAGlen=zeros(nLS,1);
        IBarrierC=zeros(nLS,1) ; IBarrierAGlen=zeros(nLS,1);
        parfor JJ=1:nLS
            
            Cplot=kk_proj(C0+gamma_StepVector(JJ)*dJdCsearch,upC,lowC);
            AGlenplot=kk_proj(AGlen0+gamma*dJdAGlensearch,upA,lowA);
            
            [Jplot(JJ),Idata(JJ),IRegC(JJ),IRegAGlen(JJ),IBarrierC(JJ),IBarrierAGlen(JJ)]=...
                CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlenplot,Cplot,n,m,alpha,rho,rhow,g,GF,Priors,Meas);
        end
        
        fVector(ig+1:ig+nLS,1)=Jplot;
        fVector(ig+1:ig+nLS,2)=Idata ;
        fVector(ig+1:ig+nLS,3)=IRegC;
        fVector(ig+1:ig+nLS,4)=IRegAGlen;
        fVector(ig+1:ig+nLS,5)=IBarrierC;
        fVector(ig+1:ig+nLS,6)=IBarrierAGlen;
        gamma_Vector(ig+1:ig+nLS)= gamma_StepVector;
        
        ind=find(~isnan(gamma_Vector)) ; gamma_Vector=gamma_Vector(ind) ;
        
        temp=fVector;
        fVector=zeros(length(ind),6);
        fVector(:,1)=temp(ind,1); %  J cost function
        fVector(:,2)=temp(ind,2); %  Idata
        fVector(:,3)=temp(ind,3); %  IRegC
        fVector(:,4)=temp(ind,4); %  IRegAGlen
        fVector(:,5)=temp(ind,5); %  IBarrierC
        fVector(:,6)=temp(ind,6); %  IBarrierAGlen
        
        [gamma_Vector,ind]=sort(gamma_Vector) ;
        fVector(:,1)=fVector(ind,1); fVector(:,2)=fVector(ind,2);
        fVector(:,3)=fVector(ind,3);fVector(:,4)=fVector(ind,4);
        fVector(:,5)=fVector(ind,5); fVector(:,6)=fVector(ind,6);
        
   
        
        figure; hold off ;
        plot(gamma_Vector,fVector(:,1),'-ob','LineWidth',2)
        hold on
        plot(gamma_Vector,fVector(:,2),'-xg')
        plot(gamma_Vector,fVector(:,3),'-xr')
        plot(gamma_Vector,fVector(:,4),'-xc')
        plot(gamma_Vector,fVector(:,5),'-xy')
        plot(gamma_Vector,fVector(:,6),'-xm')
        legend('J','DataMisfit','SystemNormC','SystemNormA','BarrierC','BarrierAGlen')
        hold on ; plot(gammaAdjoint,JMin,'*r') ; title(sprintf(' Iteration %-i ' ,iteration)) ; xlabel('\gamma') ; ylabel(' J ')
        gamma_eps=gammaAdjoint/5;
        plot([0,gamma_eps],[J0,J0-Slope0*gamma_eps],'r','LineWidth',2)
        
        SlopeAGlenBarrier0=-dIdAGlenbarrier'*dJdAGlensearch;
        plot([gamma_Vector(2),gamma_eps],[fVector(2,6),fVector(2,6)-SlopeAGlenBarrier0*gamma_eps],'r','LineWidth',2)
        
        SlopeAGlenReg=-dIdAGlenreg'*dJdAGlensearch;
        plot([gamma_Vector(2),gamma_eps],[fVector(2,4),fVector(2,4)-SlopeAGlenReg*gamma_eps],'r','LineWidth',2)
        
        if isAgrad
            SlopeAGlendata=-dIdAGlendata'*dJdAGlensearch;
            plot([gamma_Vector(2),gamma_eps],[fVector(2,2),fVector(2,2)-SlopeAGlendata*gamma_eps],'r','LineWidth',2)
        end
        
        SlopeCBarrier0=-dIdCbarrier'*dJdCsearch;
        plot([gamma_Vector(2),gamma_eps],[fVector(2,5),fVector(2,5)-SlopeCBarrier0*gamma_eps],'r','LineWidth',2)
        
        SlopeCReg=-dIdCreg'*dJdCsearch;
        plot([gamma_Vector(2),gamma_eps],[fVector(2,3),fVector(2,3)-SlopeCReg*gamma_eps],'r','LineWidth',2)
        
        if isCgrad
            SlopeCdata=-dIdCdata'*dJdCsearch;
            plot([gamma_Vector(2),gamma_eps],[fVector(2,2),fVector(2,2)-SlopeCdata*gamma_eps],'r','LineWidth',2)
        end
        
        
        
    end
%%   
    fprintf('\n ++++ Inv. it. (end): %-i \t J/J0=%-g  \t J=%-g  \t J0=%-g \t Idata=%-g \t IRegC=%-g \t IBarrierC=%-g \t IRegAGlen=%-g \t IBarrierAGlen=%-g  \n \n ',...
        iteration+Info.InverseIterations,JMin/J0,JMin,J0,IdataMin,IRegCmin,IBarrierCmin,IRegAGlenmin,IBarrierAGlenmin)
    iJ=iJ+1;  Info.JoptVector(iJ,1)=JMin; Info.JoptVector(iJ,2)=IdataMin;
    Info.JoptVector(iJ,3)=IRegCmin; Info.JoptVector(iJ,4)=IRegAGlenmin;
    Info.JoptVector(iJ,5)=IBarrierCmin; Info.JoptVector(iJ,6)=IBarrierAGlenmin;
    Info.JoptVector(iJ,7)=gammaAdjoint;
    Info.AdjointGrad{iJ}=CtrlVar.AdjointGrad;
    %%
    
    %Cest=Ctest;
    %AGlenEst=AGlentest;
    Cest=kk_proj(C0+gammaAdjoint*dJdCsearch,upC,lowC);
    AGlenEst=kk_proj(AGlen0+gammaAdjoint*dJdAGlensearch,upA,lowA);
    
    fprintf(' \t  \t \t  a=%-g \t  b=%-g \t  c=%-g  \t gamma_Min=%-g \t gamma_MinEstimate=%-g \n',gamma_a,gamma_b,gamma_c,gammaAdjoint,gamma_MinEstimate)
    
    
    %%
    
    if CtrlVar.InfoLevelAdjoint > 10000
        
        filename=sprintf('%s-AdjointIteration%04i',CtrlVar.Experiment,iteration);
        fprintf('Saving current estimate of C in %s \n',filename)
        save(filename,'CtrlVar','MUA','Cest','dJdC')
        
    end
    
    if gammaAdjoint==0 ;
        fprintf(' gamma returned equal to zero. line search has stagnated. breaking out \n')
        break
    end
    
end

I=isnan(Info.JoptVector(:,1)) ; Info.JoptVector(I,:)=[];
Info.InverseIterations=Info.InverseIterations+iteration;

end

