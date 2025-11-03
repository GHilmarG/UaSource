function [Cest,AGlenEst,JoptVector,ub,vb,ud,vd,lx,ly,gammaAdjoint]=GHGprojbfgsAdjoint(...
    CtrlVar,MUA,JoptVector,...
    s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,...
    sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
    AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
    Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
    GF)

% [C,ub,vb,JoptVector] = GHGprojbfgsAdjoint(CtrlVar,...
%         s,S,B,h,ub,vb,uMeas,vMeas,wMeas,Cd,CC,CAGlen,C_prior,AGlen_prior,coordinates,connectivity,Boundary,nip,AGlen,C,...
%         Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,n,m,alpha,rho,rhow,g,GF,DTxy,TRIxy,JoptVector)


    % Projected BFGS update
    
    %
    %  It is generally a good idea to scale the problem so that the order-of-magnitude for the typical step size is about 1.
    %  To do this do a few iterations with xscale=1 and if typical step size is alpha then
    %  set xscale=sqrt(alpha). There seems to be little need to scale the misfit funciton itself, i.e. just set fscale=1;
    %
    %  J(C)=J(x/xScale) -> C=xScale * x   and x=C/xScale
    %  dJ/dC= dJ/dx   * dx/dC  = dJ/dx   / xScale  
    %
    %  Quadradic model:  J(C)= J(C0)+  dJ/dC  (C-C0)  + (C-C0)  H (C-C0)/2
    %
    %   J(C)= J(C0)+  dJ/dC  (C-C0)  + (C-C0)  H (C-C0)/2
    %   
    %   line search:  J(C0+lambda dJ/dC) 
    %
    %
    %
    
    % This is modification of a code by C. T Kelley
    % -back tracking uses polynomial fit with an inital step estimate based on gradient
    % -vectorisation of a few loops
    % -estimates for H0 and scalings for function (fScale) and function argument (xScale)
    
    %
    % C. T. Kelley, June 11, 1998
    %
    % This code comes with no guarantee or warranty of any kind.
    %
    
    %
    % projected BFGS with Armijo rule, simple linesearch
    %
    %
    % Input: x0 = initial iterate
    %        f = objective function,
    %            the calling sequence for f should be
    %            [fout,gout]=f(x) where fout=f(x) is a scalar
    %              and gout = grad f(x) is a COLUMN vector
    %        up = vector of upper bounds
    %        low = vector of lower bounds
    %        tol = termination criterion norm(grad) < tol
    %              optional, default = 1.d-6
    %        maxit = maximum iterations (optional) default = 1000
    
    
    
    %  quadradic model
    %  fc=fc0+gc0'*(x-x0)+(x-x0)'*ddfxx*(x-x0)/2
    %  dfdx(x0)=gc0
    %  dfdx=gc0+ddfxx*(x-x0)   =>  dfdx=0 when gc=-ddfxx*(x-x0) ie dx=-ddfxx\gc
    %
    %  fc=kappa*fc0 -> kappa fc0= fc0+gc0'*gamma*dx , with dx=gc0 (i.e. steepest gradient with gamma as a step size)
    %             -> kappa fc0 =fc0 +gc0'*gc0*gamma
    %             -> gamma= (1-kappa)*fc0 /gc0'*gc0
    %
    %  if I go in the direction dsd then dx=gamma*dsd
    %  fc=fc0+gc0'*gamma*dsd
    %
    %  |fc-fc0|/gamma=gc0'*dsd  % so gc0'*dsd  is the slope with respect to gamma
    %
    %
    %
    
    
    
    nIt=CtrlVar.MaxAdjointIterations;
    if CtrlVar.AdjointRestart==0;
        JoptVector=zeros(nIt+1,6)+NaN; iJ=0;
    else
        iJ=size(JoptVector,1)-1;
        JoptVector=[JoptVector;zeros(nIt,6)+NaN];
        if iJ==-1 ; JoptVector=zeros(nIt+1,6)+NaN; iJ=0; end
    end
    
    
    
    
    iA=strfind(CtrlVar.AdjointGrad,'A'); iC=strfind(CtrlVar.AdjointGrad,'C');
    CtrlVar.isAgrad=~isempty(iA); CtrlVar.isCgrad=~isempty(iC);
    
    up=zeros(length(C),1)+CtrlVar.Cmax;
    low=zeros(length(C),1)+CtrlVar.Cmin;
    ndim=length(C);
    
    
     Itime=0;
    
    
    %% Define additional input parameters
    maxit=CtrlVar.MaxAdjointIterations;
    CtrlVar.SolveForwardModelInAdjointGradient=0;
    
    
    lsReduction=1.d0; %  controls required reducion in line search
    fRelError=1e-3;  % controls if new function values are significantly different from others
    nsmax=CtrlVar.Maximum_Number_of_BFGS_updates;  %  number of vector iterations included in BFGS update
    
    %lambdaOptEstimate=CtrlVar.AdjointInitialStep;
    
    
    %fScale=CtrlVar.AdjointfScale; % the scaling of the cost function is now done in CostFunctionValueAndGradient
                                  % but fscale is still needed for scaling the gradient
    fScale=1;
    xScale=CtrlVar.AdjointxScale; H0=zeros(ndim,1)+1;
    %xScale=1 ; H0=zeros(ndim,1)+CtrlVar.AdjointxScale; % estimate of initial inverse Hessian
    
    
    fVector=zeros(100,1)+NaN ;  lambdaVector=zeros(100,1)+NaN ; % data needed for plotting purposes
    
    %% Scale input variables
    % if f(x) is the user supplied function, I minimize g(x)=fscale f(xScale x),
    % so grad g = fScale xScale grad f; g(x0+alpha grad g)= fscale f(xScale x0 + fscale xScale^2 grad f)
    %H0=H0/(fScale*xScale*xScale);
    x0=C/xScale;  up=up/xScale ; low=low/xScale;
    
    %%
    
    xc=x0;
    kku=up ; kkl=low;
    if any(kkl>kku) ; error(' lower bound exceeds upper bound ') ; end  % GHG
    
    % put initial iterate in feasible set
    %
    if norm(xc - kk_proj(xc,kku,kkl)) > 0
        disp(' initial iterate not feasible ');
        xc=kk_proj(xc,kku,kkl);
    end
    
    if CtrlVar.AdjointRestart==1
        load BFGSdata ystore sstore itc numf numg numh ithist ns
        fprintf(' Loaded BFGS from a previous run because this is an adjoint restart run \n')
        maxit=CtrlVar.MaxAdjointIterations;
        ithist=[ithist;zeros(maxit,6)];
    else
        ystore=zeros(ndim,nsmax); sstore=ystore; ns=0;
        itc=1; numf=1; numg=1; numh=0;
        ithist=zeros(maxit,6);
    end
    
    
    C=xScale*xc;
    
%     [fc,gc,Idata,IReg,IBarrier] = CostFunctionValueAndGradient(s,S,B,h,ub,vb,coordinates,connectivity,Boundary,nip,...
%         AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
%         LAdjoint,LAdjointrhs,lambdaAdjoint,DTxy,TRIxy,...
%         uMeas,vMeas,wMeas,C_prior,AGlen_prior,Cd,CAGlen,CC,GF,CtrlVar);
%     
[fc,gc,Idata,IReg,IBarrier] = CostFunctionValueAndGradient(CtrlVar,MUA,JoptVector,...
    s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,...
    sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
    AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
    Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
    GF);


    
    iJ=iJ+1; JoptVector(iJ,:)=[fc,Idata,IReg,IBarrier,norm(gc),0];
    
    gc=fScale*xScale*gc;
    
    kappa=0.5 ; lambdaOptEstimate=(1-kappa)*fc/((H0.*gc)'*(H0.*gc));
    %lambdaOptEstimate=1;
    fprintf(' initial lambda opt estimate %-g \n ',lambdaOptEstimate)
    fprintf(' initial mistfit=%-g \t with slope %-g \n ',fc,sqrt(gc'*gc))
    
    
    
    %pgc=xc - kk_proj(xc - gc,kku,kkl);
    pgc=xc - kk_proj(xc - lambdaOptEstimate*H0.*gc,kku,kkl);  % GHG, this is the projected step in direction gc
                                                              % gc is the steepest gradient
    
    
    %ia=0; alist=zeros(ndim,1);
    tst=kku-kkl;
    lim1=.5*min(tst);
    
    epsilon=min(lim1,norm(pgc));
    fprintf(' epsilon set at %-g with lim1=%-g and norm(pgc)=%-g \n ',epsilon,lim1,norm(pgc))
    ia=sum(xc==kku | xc==kkl); % GHG,  number of active constraints
    
    alist=min(kku-xc,xc-kkl)< epsilon ; % GHG
    fprintf(' initial fraction of epsilon active constrains %-g \n ',sum(alist)/length(alist))
    
    
    ithist(itc,1) = fc;
    ithist(itc,2) = norm(pgc);
    ithist(itc,3)=0;
    ithist(itc,4)=itc-1;
    ithist(itc,5)=ia/ndim;
    ithist(itc,6)=-gc'*gc;
    
    %%
    iteration=0; iytsnegative=0;
    while( iteration < maxit)
        iteration=iteration+1;
        
        
        fprintf('\n ##############     itc=%-i \t ns=%-i  \t iteration=%-i \t maxitc=%-i  ############  \n ',itc,ns,iteration,maxit)
        ifVector=1;
        fVector(ifVector)=fc; lambdaVector(ifVector)= 0 ; ifVector=ifVector+1;
        xc0=xc; fc0=fc; gc0=gc; % function value and gradient at the beginning of step
        
        
        %dsd=-gc;
        dsd=-H0.*gc;   %  Newton gives dsd=-H0.*gc,  here  H0 is an approximation of the inverse Hessian 
        dsd=bfgsrp(ystore,sstore,ns,dsd,alist,H0); % dsd(i)=0 if alist(i)=1
        
        % dsd=dsd+proja(-gc,alist); % dsd(i)=-gc(i) if alist(i)=1, dsd(i) given by bfgsrp otherwise
        dsd=dsd+proja(-H0.*gc,alist); % dsd(i)=-gc(i) if alist(i)=1, dsd(i) given by bfgsrp otherwise
        % if all constraints are estimated to be active then dsd=-H0.*gc (test
        % this)
        
        %lambdaOptEstimate=lambdaMin;
        kappa=0.5 ;   %  gc the gradient, dsd the search direction
        
        lambdaOptEstimate=(1-kappa)*fc/(-gc'*dsd);
        lambda=lambdaOptEstimate;
        %lambdaOptEstimate=1;
        fprintf('  lambda opt estimate %-g \n ',lambdaOptEstimate)
        
        
        xt=kk_proj(xc+lambda*dsd,kku,kkl); % using full Newton step takes us to xt (stepsize: lambda, direction: dsd)
        C=xScale*xt;
        

        ft = CostFunctionValueAndGradient(s,S,B,h,ub,vb,coordinates,connectivity,Boundary,nip,...
            AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
            LAdjoint,LAdjointrhs,lambdaAdjoint,DTxy,TRIxy,...
            uMeas,vMeas,wMeas,C_prior,AGlen_prior,Cd,CAGlen,CC,GF,CtrlVar);
        
        
        
        fVector(ifVector)=ft; lambdaVector(ifVector)= lambda ; ifVector=ifVector+1;
        %ft=feval(f,xt);
        numf=numf+1;
        
        
        iarm=0; itc=itc+1;
        pl=xc - xt; % Newton step size
        fgoal=fc-(gc'*pl)*lsReduction; % if ft<=fgoal, the full step is accepted and no line search performed
        
        fprintf(' Newton step %-i \t  fc=%-g  \t ft=%-g \t fgoal=%-g \t ft new/ft last=%-g \t ft new/fc=%-g \t lambda=%-g \n ',iarm,fc,ft,fgoal,ft/fc,ft/fc,lambda)
        
        if ft/fc < 0.5 ;
            fprintf(' fc=%-g \t ft=%-g \t fgoal=%-g and full step is accepted with ft/fc=%-g \n ',fc,ft,fgoal,ft/fc)
            NewtonAccepted=1;
            fc=ft; lambdaMin=lambda;
        else
            fprintf(' fc=%-g \t ft=%-g \t fgoal=%-g and full step not accepted \n ',fc,ft,fgoal)
            NewtonAccepted=0;
        end
        
        
        
        
        %         fprintf(' exploratory calculations  ')
        %         tStart=tic;
        %         ifVector=ifVector-1;
        %         lambdaPlot=logspace(-3,3,8);
        %         parfor I=1:8
        %
        %             CPlot=xScale*kk_proj(xc0+lambdaPlot(I)*dsd,kku,kkl);
        %             [uPlot,vPlot]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,CPlot,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
        %             fPlot=MisfitFunction(uPlot,vPlot,wint,uMeas,vMeas,wMeasInt,CPlot,C_prior,AGlen,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar); fPlot=fScale*fPlot;
        %             fVector(ifVector+I)=fPlot; lambdaVector(ifVector+I)= lambdaPlot(I) ;
        %
        %         end
        %         ifVector=ifVector+8+1;
        %         tElapsed=toc(tStart); fprintf(' done in %-f sec \n',tElapsed)
        %
        %% line search
        q0=fc; qp0=gc'*dsd ; bhigh=.5; blow=.1; exFactor=5;
        
        
        
        %% extrapolation step
        % if ft/fc>0.25 and I have less than 10 updates to the Hessian do extrapolation until past minimum
        if NewtonAccepted==0 % line search
            
            iExStep=0; iExStepMax=7; lamc=lambda; qc=ft; qm=ft; qmm=qm;
            while iExStep<=iExStepMax && ft/fc>0.5
                
                xtLast=xt;
                iExStep=iExStep+1; iarm=iarm+1;
                
                
                lambda=exFactor*lambda;
                
                xt=kk_proj(xc+lambda*dsd,kku,kkl);
                C=xScale*xt;
                
                %                 [u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
                %                 ft=MisfitFunction(u,v,wint,uMeas,vMeas,wMeasInt,C,C_prior,AGlen,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar); ft=fScale*ft;
                %                  ft=fScale*ft;
                %
                ft = CostFunctionValueAndGradient(s,S,B,h,ub,vb,coordinates,connectivity,Boundary,nip,...
                    AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
                    LAdjoint,LAdjointrhs,lambdaAdjoint,DTxy,TRIxy,...
                    uMeas,vMeas,wMeas,C_prior,AGlen_prior,Cd,CAGlen,CC,GF,CtrlVar);
                
                
                
                qm=qc; lamm=lamc; lamc=lambda; qc=ft;
                numf = numf+1;
                
                fVector(ifVector)=ft; lambdaVector(ifVector)= lambda ; ifVector=ifVector+1;
                fprintf(' Extrapolation step %-i \t  fc=%-g  \t ft=%-g \t fgoal=%-g \t ft new/ft last=%-g \t ft new/fc=%-g \t lambda=%-g \n ',iExStep,fc,ft,fgoal,qc/qm,qc/fc,lambda)
                
                
                if ft> qm  %  if past minimum  exit, but first do cubic/parbolic fit
                    ft=qm; xt=xtLast; lambda=lamm;
                    
                    lambdaTest=NaN;
                    if ifVector>4;
                        lambdaTest= parabolamin(lambdaVector(ifVector-3:ifVector-1),fVector(ifVector-3:ifVector-1) ) ;
                        if isnan(lambdaTest)
                            lambdaTest=max([blow,lambdaTest]);
                            lambdaTest=min([exFactor*max(lambdaVector),lambdaTest]);
                        end
                        fprintf(' parabolic fit through last three points in extrapolation step gives lambda=%-g \n ',lambdaTest)
                    end
                    
                    if isnan(lambdaTest)
                        lambdaTest=GHGpolymod(q0, qp0, lamc, qc, blow, bhigh, lamm, qm);
                        fprintf(' cubic fit trough last three points and slope at zero in exprapolation step gives lambda=%-g \n ',lambdaTest)
                    end
                    
                    xtTest=kk_proj(xc+lambdaTest*dsd,kku,kkl);
                    C=xScale*xtTest;
                    if ~isreal(C) ; save TestSave C lambdaTest dsd xc kku kkl q0 qp0 lamc qc lamm qm ; error(' C not real') ;  end
                    
                    %[u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
                    %ftTest=MisfitFunction(u,v,wint,uMeas,vMeas,wMeasInt,C,C_prior,AGlen,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);  ftTest=fScale*ftTest;
                    
                    ftTest = CostFunctionValueAndGradient(s,S,B,h,ub,vb,coordinates,connectivity,Boundary,nip,...
                        AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
                        LAdjoint,LAdjointrhs,lambdaAdjoint,DTxy,TRIxy,...
                        uMeas,vMeas,wMeas,C_prior,AGlen_prior,Cd,CAGlen,CC,GF,CtrlVar);
                    
                    
                    
                    qm=qc; lamm=lamc; lamc=lambdaTest; qc=ftTest;
                    fVector(ifVector)=ftTest; lambdaVector(ifVector)= lambdaTest ; ifVector=ifVector+1;
                    if ftTest<ft % test if value is indeed smaller and accept if so
                        lambda=lambdaTest ; ft=ftTest ; xt=xtTest;
                        fprintf(' accepting ftTest in extrapolation step with ft=%-g \t ft/fc=%-g \t lambda=%-g \n ',ft,ft/fc,lambda)
                    end
                    break
                end
            end % extrapolation step
            %%
            
            %% backtracking step
            
            while (ft > fgoal || (iarm==0  && ft/fc>0.5 )) && qc < qmm % backtracking as long as: 1) ft>ftgoal or new value smaller than last one
                iarm=iarm+1;
                
                
                if iarm==1
                    bhighQuad=2.5; % start with a general quadradic fit, ie not forcing backtracking
                    lambda=GHGpolymod(q0, qp0, lamc, qc, blow, bhighQuad);
                else
                    lambda=GHGpolymod(q0, qp0, lamc, qc, blow, bhigh, lamm, qm);
                end
                
                xt=kk_proj(xc+lambda*dsd,kku,kkl);
                pl=xc-xt;
                
                C=xScale*xt;
                
                % [u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
                % ft=MisfitFunction(u,v,wint,uMeas,vMeas,wMeasInt,C,C_prior,AGlen,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
                
                
                
                ft = CostFunctionValueAndGradient(s,S,B,h,ub,vb,coordinates,connectivity,Boundary,nip,...
                    AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
                    LAdjoint,LAdjointrhs,lambdaAdjoint,DTxy,TRIxy,...
                    uMeas,vMeas,wMeas,C_prior,AGlen_prior,Cd,CAGlen,CC,GF,CtrlVar);
                
                
                
                %ft=feval(f,xt);
                numf = numf+1;
                
                qmm=qm ; qm=qc; lamm=lamc; lamc=lambda; qc=ft;
                
                
                fVector(ifVector)=ft; lambdaVector(ifVector)= lambda ; ifVector=ifVector+1;
                
                if(iarm > iExStep+20)
                    disp(' Armijo error in steepest descent \n ')
                    break;
                elseif (iarm > iExStep+5) && abs(1-ft/fc)<fRelError && abs(1-qm/qc)<fRelError
                    % allowing for a slight increase in cost function if I have
                    % already done three backtrack steps solution if very close to start
                    %
                    break
                end
                fgoal=fc-(gc'*pl)*lsReduction;
                fprintf(' Backtracking step %-i \t  fc=%-g  \t ft=%-g \t fgoal=%-g \t ft new/ft last=%-g \t ft new/fc=%-g \t lambda=%-g \n ',iarm,fc,ft,fgoal,qc/qm,qc/fc,lambda)
            end
            
            
            [~,ind]=min(fVector(2:end)) ; lambdaMin=lambdaVector(ind+1);
            xt=kk_proj(xc+lambdaMin*dsd,kku,kkl); C=xScale*xt;
            
            
            numf=numf+1;
            
        end
        
        [fc,gp,Idata,IReg,IBarrier] = CostFunctionValueAndGradient(s,S,B,h,ub,vb,coordinates,connectivity,Boundary,nip,...
            AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
            LAdjoint,LAdjointrhs,lambdaAdjoint,DTxy,TRIxy,...
            uMeas,vMeas,wMeas,C_prior,AGlen_prior,Cd,CAGlen,CC,GF,CtrlVar);
        
        iJ=iJ+1; JoptVector(iJ,:)=[fc,Idata,IReg,IBarrier,norm(gp),0];
        
        gp=fScale*xScale*gp;
        
        
        %figure ; trisurf(TRIxy,coordinates(:,1),coordinates(:,2),gc) ;xlabel('x') ; ylabel('y') ; title(sprintf('gradient at itc=%-i ',itc))
        
        
        %[gp]=feval(fgrad,xt);
        numg=numg+1;
        
        yLocal=gp-gc; % difference of the gradient vectors at new and last starting point
        sLocal=xt-xc; %
        
        gc=gp; xc=xt;
        %pgc=xc-kk_proj(xc-gc,kku,kkl);  % this seems to be size of next projected full Newton step with Hess=1
        
        pgc=xc-kk_proj(xc-lambdaMin*H0.*gc,kku,kkl);  % GHG
        epsilon=min(lim1,norm(pgc));
        %epsilon=1e-3;
        alist=min(kku-xc,xc-kkl)< epsilon  ;
        fprintf(' fraction of epsilon active constrains %-g \n ',sum(alist)/length(alist))
        
        
        alphabb=yLocal'*yLocal/(sLocal'*yLocal);   fprintf(' alphabb=%-g \n ',alphabb)
        
        yLocal=proji(yLocal,alist); sLocal=proji(sLocal,alist);
        
        %
        %   restart if y'*s is not positive or we're out of room
        %
        % y=gp-gc is grad f(xt)- grad f(xc) where xc starting point
        % s=xt-xc=xc+lambdaMin*dsd-xc where xc starting point
        % i.e. (grad f(xt)-grad f(xc))' * (xt-xc)
        
        
        yts=yLocal'*sLocal;  % positive curvature condition
        
        if ns>=nsmax
            ns=0;
            fprintf(' resetting BFGS update because maximum of iterations reached')
        elseif yts <=0 ;
            iytsnegative=iytsnegative+1;
            ns=0;
            fprintf(' yts=%-g <=0 , resetting BFGS update \n ',yts)
        elseif yts > 0
            iytsnegative=0;
            ns=ns+1;
            alphaLocal=sqrt(yts);
            sstore(:,ns)=sLocal/alphaLocal; ystore(:,ns)=yLocal/alphaLocal;
        end
        ithist(itc,2)=norm(pgc); ithist(itc,1) = fc;  ithist(itc,4)=itc-1; ithist(itc,3)=iarm;
        %ia=0; for i=1:ndim; if(xc(i)==kku(i) || xc(i)==kkl(i)) ; ia=ia+1; end; end;
        ia=sum(xc==kku | xc==kkl); % GHG
        
        ithist(itc,5)=ia/ndim; ithist(itc,6)=-gc'*dsd;
        
        
        %% trust region
        % update the upper estimate for lambda
        % is actual reduction similar to the one expeced from quadradic model?
        
        ratiof=2*(fc-q0)/(gc0'*dsd);
        ratioLambda=lambdaMin/lambdaOptEstimate;
        fprintf(' reduction=%-g \t ratio between actual and expected reduction in f for a full Newton step %-g \n ',...
            fc/fc0,ratiof)
        fprintf(' lambdaMin=%-g \t  ratio between optimal stepsize and upper estimate for stepsize is %-g \n ',...
            lambdaMin,ratioLambda)
        
        
        
        
        if CtrlVar.InfoLevelAdjoint>=100;
            fprintf(' calculations needed for plotting purposes ')
            tStart=tic;
            ifVector=ifVector-1;
            parfor I=1:8
                
                lambdaPlot=lambdaMin/1000+(I-2)*lambdaMin/4;
                CPlot=xScale*kk_proj(xc0+lambdaPlot*dsd,kku,kkl);
                
                %[uPlot,vPlot]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,CPlot,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
                %fPlot=MisfitFunction(uPlot,vPlot,wint,uMeas,vMeas,wMeasInt,CPlot,C_prior,AGlen,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);fPlot=fScale*fPlot;
                
                
                fPlot = CostFunctionValueAndGradient(s,S,B,h,ub,vb,coordinates,connectivity,Boundary,nip,...
                    AGlen,CPlot,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,...
                    LAdjoint,LAdjointrhs,lambdaAdjoint,DTxy,TRIxy,...
                    uMeas,vMeas,wMeas,C_prior,AGlen_prior,Cd,CAGlen,CC,GF,CtrlVar);
                
                
                fVector(ifVector+I)=fPlot; lambdaVector(ifVector+I)= lambdaPlot ;
                
            end
            
            tElapsed=toc(tStart); fprintf(' done in %-f sec \n',tElapsed)
            
            lambdaVector=lambdaVector(~isnan(lambdaVector)) ; fVector=fVector(~isnan(lambdaVector));
            [lambdaVector,ind]=unique(lambdaVector) ; fVector=fVector(ind);
            [lambdaVector,ind]=sort(lambdaVector) ; fVector=fVector(ind);
            
            
            
            if CtrlVar.doplots==1
                figure ; plot(lambdaVector,fVector,'-+') ; hold on ; plot(lambdaMin,fc,'go')
                plot([lambdaVector(2),lambdaVector(3)],[fVector(2),fVector(2)+(lambdaVector(3)-lambdaVector(2))*qp0],'-rs')
                title(sprintf('J at itc=%-i for ns=%0i ',itc-1,ns)) ; xlabel('\lambda')
            end
            
            sol=[ones(4,1) lambdaVector(2:5) lambdaVector(2:5).*lambdaVector(2:5)]\fVector(2:5);
            fprintf(' Slope calculated=%-g \t slope estimated=%-g \t %-g \n',...
                sol(2)+2*sol(3)*lambdaVector(1),gc0'*dsd)
            
            fVector=zeros(100,1)+NaN ;  lambdaVector=zeros(100,1)+NaN ;
        end
        
        
    end
    
    
    
    
    x=xc*xScale ;
    C=x;
    
    [ub,vb]=SSTREAM2dNR(s,S,B,h,ub,vb,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
    
    save BFGSdata ystore sstore itc numf numg numh ithist ns
    
end


function px = kk_proj(x,kku,kkl)
    %
    % projection onto active set
    %
    %ndim=length(x);
    %px=zeros(ndim,1);
    px=min(kku,x);
    px=max(kkl,px);
    
end

function dnewt=bfgsrp(ystore,sstore,ns,dsd,alist,H0)
    %
    % bfgsrp
    %
    % C. T. Kelley, Dec 20, 1996
    %
    % This code comes with no guarantee or warranty of any kind.
    %
    % This code is used in projbfgs.m
    %
    % There is no reason to ever call this directly.
    %
    % form the product of the generalized inverse of the
    % bfgs approximate Hessian
    % with a vector using the recursive approach
    %
    dnewt=proji(dsd,alist);
    if (ns==0)
        H0=projhess(H0,alist); dnewt=H0.*dnewt;
        fprintf(' \n ')
        return
    end
    sstore(:,ns)=proji(sstore(:,ns),alist);
    ystore(:,ns)=proji(ystore(:,ns),alist);
    beta=sstore(:,ns)'*dsd; dnewt=dsd-beta*ystore(:,ns);
    ndim=length(dsd); xlist=zeros(ndim,1);
    dnewt=bfgsrp(ystore,sstore,ns-1,dnewt,xlist,H0);
    dnewt=dnewt+(beta-ystore(:,ns)'*dnewt)*sstore(:,ns);
    dnewt=proji(dnewt,alist);
    fprintf(' %-i ',ns)
    
end

function px=proji(x,alist)
    %
    % projection onto epsilon-inactive set
    %
    % sets elements of x to zero where corresponding elements of alist iare equal to on
    %
    %ndim=length(x);
    px=x;
    %for k=1:ndim ; if alist(k) == 1 ; px(k)=0; end; end
    
    px(alist==1)=0; % GHG
    
end


function px=proja(x,alist)
    %
    %  projection onto epsilon-active set
    %
    %ndim=length(x);
    px=x;
    %for k=1:ndim ; if alist(k) == 0 ; px(k)=0; end; end
    px(alist==0)=0; % GHG
    
end


function px=projhess(x,alist)
    %
    % projection onto epsilon-inactive set
    %
    %ndim=length(x);
    px=x;
    %for k=1:ndim ; if alist(k) == 1 ; px(k)=0; end; end
    
    px(alist==1)=0; % GHG
    
end
