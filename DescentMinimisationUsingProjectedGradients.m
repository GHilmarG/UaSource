
function [C,AGlen,u,v,JoptVector]=DescentMinimisationUsingProjectedGradients(sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
        s,S,B,h,u,v,coordinates,connectivity,Xint,Yint,xint,yint,Boundary,DTxy,TRIxy,DTint,TRIint,...
        nip,AGlen,C,Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
        n,m,alpha,rho,rhow,g,CtrlVar,Itime,JoptVector)
    
    CtrlVar.SolveForwardModelInAdjointGradient=1;
    x=coordinates(:,1); y=coordinates(:,2);
    dJdC=0; dJdescentC=0;   dJdAGlen=0; dJdescentA=0;
    
   
    
    wMeasInt=Grid1toGrid2(DTxy,wMeas,Xint,Yint);
    
    for iteration=ItStart:ItStart-1+CtrlVar.MaxAdjointIterations
        close all
        
        
        xvector=zeros(100,1)+NaN; fvector=zeros(100,1)+NaN;
        
        %C(C<1000*CtrlVar.Cmin)=1000*CtrlVar.Cmin;
        
        % get gradient of misfit with respect to C
        
        
        dJdescentClast=dJdescentC; dJdescentAlast=dJdescentA;
        
        
        [dJdescentC,dJdescentA,u,v,wint,lx,ly]=AdjointGradientNR2d(s,S,B,h,u,v,uMeas,vMeas,wMeasInt,Cd,CC,CAGlen,C_prior,AGlen_prior,coordinates,connectivity,Boundary,nip,AGlen,C,...
            Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,n,m,alpha,rho,rhow,g,Itime,CtrlVar,iteration,DTxy,TRIxy);
        
        
        
        
        if iteration>(2+ItStart) && CtrlVar.AdjointConjugatedGradients==1;
            
            fprintf(' conj grad step \n')
            
            [dJdC,ConjGradAngleC]=NewConjugatedGrad(dJdescentC,dJdescentClast,dJdC);
            [dJdAGlen,ConjGradAngleA]=NewConjugatedGrad(dJdescentA,dJdescentAlast,dJdAGlen);
            
            fprintf(' ConjGradAngleC=%-g degrees, \t ConjGradAngleA=%-g degrees \n',ConjGradAngleC,ConjGradAngleA)
            
            
        else
            dJdC=-dJdescentC; ConjGradAngleC=0; dJdAGlen=-dJdescentA; ConjGradAngleA=0;
        end
        
        
        
        
        
        % f=f(x0)+f(x0+gamma* grad f) hence df/dgamma=grad f * grad f
        %
        %
        %
        SlopeGamma=-dJdC'*dJdC;
        
        figure; trisurf(TRIxy,x,y,lx) ;  title(sprintf(' lx at iteration %-i ',iteration) )
        figure; trisurf(TRIxy,x,y,ly) ;  title(sprintf(' ly at iteration %-i ',iteration) )
        
        figure; trisurf(TRIxy,x,y,dJdC) ;  title(sprintf(' dJdC at iteration %-i ',iteration) )
        figure ; trisurf(TRIxy,x,y,dJdAGlen) ;  title(sprintf(' dJdAGlen at iteration %-i ',iteration) )
        figure ; trisurf(TRIxy,x,y,u) ;  title(sprintf(' u at iteration %-i ',iteration) )
        figure ; trisurf(TRIxy,x,y,v) ;  title(sprintf(' v at iteration %-i ',iteration) )
        speed=sqrt(u.*u+v.*v) ; speedMeas=sqrt(uMeas.*uMeas+vMeas.*vMeas);
        figure ; trisurf(TRIxy,x,y,speed)     ;  title(sprintf(' speed at iteration %-i ',iteration) )
        figure ; trisurf(TRIxy,x,y,speedMeas) ;  title(' measured speed  ')
        figure ; trisurf(TRIxy,x,y,speed-speedMeas) ;  title(' speed - speedMeas ')
        
        
        
        [J0,Idata0,IRegC0]=MisfitFunction(u,v,wint,uMeas,vMeas,wMeasInt,C,C_prior,AGlen,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,CtrlVar);
        
        fprintf('It #%-i, \t J0=%-g \n',iteration,J0)
        JoptVector(iteration)=J0 ;
        
        xvector(1)=0 ; fvector(1)=J0;
        
        if iteration==ItStart ;
            if CtrlVar.AdjointInitialStep>0
                gammaStep=CtrlVar.AdjointInitialStep;
            else
                gammaStep=-J0/SlopeGamma ;
            end
        end
        
        gamma=gammaStep;  % first step
        
        [J,u,v]=fJNR2d(gamma,dJdC,dJdAGlen,uMeas,vMeas,wMeasInt,s,S,B,h,u,v,wint,Cd,CAGlen,CC,AGlen_prior,C_prior,coordinates,connectivity,Boundary,nip,...
            AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
        
        
        
        
        Jmin=J0 ; ax=0 ; fa=J0 ;
        cx=gamma ; fc=J;
        xvector(2)=gamma ; fvector(2)=J;
        itLineSearch=0;
        
        % Line Search using projected gradient
        while Jmin>=J0 && itLineSearch<=2
            
            itLineSearch=itLineSearch+1;
            if CtrlVar.InfoLevelAdjoint>10; fprintf(' Line search \n') ; end
            
            fprintf(' initial bracket estimate : \n ')
            fprintf(' ax=%-g, \t \t fa=%-g \n',ax,fa)
            fprintf(' cx=%-g, \t \t fc=%-g \n',cx,fc)
            
            
            %% first make sure that minimum is bracketed between ax and cx, ie that fc=f(cx) > fa=f(ax) with cx>ax
            if fc<fa  % extrapolate to bracket minimum
                fprintf(' Extrapolate to bracket minimum \n')
                nStepIncrease=CtrlVar.AdjointExtrapolateStepFactor; nfmax=CtrlVar.AdjointMaxFuncEvalinExtrapolateStep ; ytolfactor=1.1;
                
                [cx,fc,xvectorTemp,fvectorTemp,ifunc,u,v]=LineSearchExtrapolate(@(gamma)...
                    fJNR2d(gamma,dJdC,dJdAGlen,uMeas,vMeas,wMeasInt,s,S,B,h,u,v,wint,Cd,CAGlen,CC,AGlen_prior,C_prior,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar),...
                    J0,ax,cx,nStepIncrease,nfmax,ytolfactor);
                
                
                xvector=[xvector(1);xvector(2);xvectorTemp(1:ifunc)] ; fvector=[fvector(1);fvector(2);fvectorTemp(1:ifunc)];
                
                
                [xvector,ind]=unique(xvector) ; fvector=fvector(ind);
                [xvector,ind]=sort(xvector) ; fvector=fvector(ind);
                
                for II=1:numel(xvector)
                    fprintf(' xvector=%-g,\t fvector=%-g \n',xvector(II),fvector(II))
                end
                
                [~,ind]=min(fvector) ;
                
                if ind>1 ; ax=xvector(ind-1) ; end
                if ind<numel(xvector) ; cx=xvector(ind+1); end
                
                fprintf(' ax=%-g, \t \t fa=%-g \n',ax,fa)
                fprintf(' cx=%-g, \t \t fc=%-g \n',cx,fc)
            else
                
                fprintf(' initial step to the right of minimum \n')
                fprintf(' \t ax=%-g \t cx=%-g \n',ax,cx)
                fprintf(' \t fa=%-g \t fc=%-g \n',fa,fc)
                for II=1:2
                    fprintf(' xvector=%-g,\t fvector=%-g \n',xvector(II),fvector(II))
                end
            end
            
            [Jmin,ind]=min(fvector) ; gammamin=xvector(ind);
            
            
            if Jmin>0.1*J0 % do line search if not sufficient decrease
                
                CtrlVar.fminbndTol=(cx-ax)/1e4;
                fprintf(' fminbdn \n ')
                TolFun=0.01*J0; % exit linesearch if function value less that this after five function calls
                
                [gammamin,Jmin,exitflag,fminbndOutput,u,v]=fminbndAdjoint2d(@(dummy) ...
                    fJNR2d(dummy,dJdC,dJdAGlen,uMeas,vMeas,wMeasInt,s,S,B,h,u,v,wint,Cd,CAGlen,CC,AGlen_prior,C_prior,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar),...
                    ax,cx,optimset('TolFun',TolFun,'TolX',CtrlVar.fminbndTol,'Display','final','MaxFunEvals',CtrlVar.fminbndMaxIterations));
                
                
                
                ind=~isnan(xvector) & ~isnan(fvector); xvectorTemp1=xvector(ind); fvectorTemp1=fvector(ind);
                ind=~isnan(fminbndOutput.xVector) & ~isnan(fminbndOutput.fVector); xvectorTemp2=fminbndOutput.xVector(ind); fvectorTemp2=fminbndOutput.fVector(ind);
                
                xvector=[xvectorTemp1(:); xvectorTemp2(:)]; fvector=[fvectorTemp1(:); fvectorTemp2(:)];
                
                
                
                switch exitflag
                    case 1
                        fprintf(' fminbnd: converged to a solution with gammamin=%-g and J=%-g \n',gammamin,J)
                    case 0
                        fprintf(' fminbnd: maximum number of iterations reached \n')
                    otherwise
                        fprintf(' fminbnd: exitflag %-g \n  ',exitflag)
                end
                
                fprintf(' fminbnd: \t Number of iterations %-g \t Number of function evaluations %-g \n ',fminbndOutput.iterations,fminbndOutput.funcCount)
                fprintf(' \t ax=%-g \t gammamin=%-g \t cx=%-g \n ',ax,gammamin,cx)
                fprintf(' \t fa=%-g \t Jmin=%-g \t fc=%-g \n ',fa,Jmin,fc)
                if Jmin > 1.05*J0 ;  % allow some increase
                    cx=cx/100 ;
                    fprintf(' fminbdn did not find a minumum, try a new and a smaller bracket %-g \n ',cx) ;
                    % gammamin=CtrlVar.AdjointMinStep;
                    
                    % Most likely this happened because a local minimum was
                    % found that is larger than J0, try to reduce the
                    % initial bracket
                    
                    
                    [fc,u,v]=fJNR2d(cx,dJdC,dJdAGlen,uMeas,vMeas,wMeasInt,s,S,B,h,u,v,wint,Cd,CAGlen,CC,AGlen_prior,C_prior,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
                end
            end
            
            
        end
        
        
        
        
        if CtrlVar.InfoLevelAdjoint>100;
            Nstep=10;
            ImisfitTestVector=zeros(Nstep,1)+NaN; gammaTestVector=zeros(Nstep,1)+NaN;
            
            
            for I=1:Nstep
                %gammaTest=ax+(I-1)*(cx-ax)/Nstep;
                %gammaTest=(I-1)*cx/Nstep;
                gammaTest=I*gammamin/Nstep/20;
                
                
                [Jtest,u,v]=...
                    fJNR2d(gammaTest,dJdC,dJdAGlen,uMeas,vMeas,wMeasInt,s,S,B,h,u,v,wint,Cd,CAGlen,CC,AGlen_prior,C_prior,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
                ImisfitTestVector(I+1)=Jtest;
                gammaTestVector(I+1)=gammaTest;
            end
            
            
            xvector=[xvector;gammaTestVector]; fvector=[fvector;ImisfitTestVector];
        end
        
        [xvector,ind]=unique(xvector) ; fvector=fvector(ind);
        [xvector,ind]=sort(xvector) ; fvector=fvector(ind);
        figure; plot(xvector,fvector,'-+') ;
        hold on ; plot(gammamin,Jmin,'or') ;
        if CtrlVar.InfoLevelAdjoint>100; hold on ; plot(gammaTestVector,ImisfitTestVector,'xg') ; end
        plot([xvector(1),xvector(2)],[fvector(1),fvector(2)],'r')
        plot([xvector(1),xvector(2)],[fvector(1),fvector(1)+SlopeGamma*(xvector(2)-xvector(1))],'g')
        
        hold off
        
        
        SlopeGammaNumerical=(fvector(2)-fvector(1))/(xvector(2)-xvector(1));
        
        fprintf('\n  SlopeGamma Adjoint=%-g, \t SlopeGamma Numerical=%-g \t Adjoint Slope/Numerical Slope=%-g \n \n',SlopeGamma,SlopeGammaNumerical,SlopeGamma/SlopeGammaNumerical)
        
        fprintf(' Estimating step size based on J0/SlopeGamma gives %-g \n \n ',-J0/SlopeGamma)
        
        
        % updating C and A
        
        figure ; trisurf(TRIxy,x,y,C) ;  title(sprintf('C at start of iteration %-i ',iteration))
        figure ; trisurf(TRIxy,x,y,AGlen) ;  title(sprintf('AGlen at start of iteration %-i ',iteration))
        
        Cold=C; C=ProjGradient(C,dJdC,gammamin,CtrlVar.Cmin,CtrlVar.Cmax);
        AGlenold=AGlen; AGlen=ProjGradient(AGlen,dJdAGlen,gammamin,CtrlVar.AGlenmin,CtrlVar.AGlenmax);
        
        
        %
        d=C-Cold ; d0=gammamin*dJdC;
        
        Angle=acosd(d'*d0/norm(d)/norm(d0));
        
        
        fprintf('Angle between C with and without projection is %-g \n ',Angle)
        
        
        figure ; trisurf(TRIxy,x,y,C) ;  title(sprintf('C at end of iteration %-i ',iteration))
        figure ; trisurf(TRIxy,x,y,C-Cold) ;  title(sprintf('C-Cold at iteration  %-i ',iteration))
        
        figure ; trisurf(TRIxy,x,y,AGlen) ;  title(sprintf('AGlen at end of iteration %-i ',iteration))
        figure ; trisurf(TRIxy,x,y,AGlen-AGlenold) ;  title(sprintf('AGlen-AGlenold at iteration  %-i ',iteration))
        
        
        
        
        % calculate u and v for the new values of A and C
        [u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
        
        [J1,Idata1,IRegC1]=MisfitFunction(u,v,wint,uMeas,vMeas,wMeasInt,C,C_prior,AGlen,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,CtrlVar);
        JoptVector(iteration+1)=J1 ;
        
        figure(100001) ; semilogy(0:iteration,JoptVector(1:iteration+1),'o-') ; title('J') ; xlabel('Iteration')
        gammaStep=100*gammamin; % use last step if
        
        
        fprintf('It # %-4i \t , J0=%-14.7g , J=%-14.7g , J/J0=%-14.7g, gammamin=%-14.7g, angle C=%-15.7g,  angle A=%-15.7g \n',iteration,J0,Jmin,Jmin/J0,gammamin,ConjGradAngleC,ConjGradAngleA)
        fprintf(' J0=%-g \t Idata0=%-f \t  IRegC0=%-f \n ',J0,Idata0,IRegC0)
        fprintf(' J1=%-g \t Idata1=%-f \t  IRegC1=%-f \n ',J1,Idata1,IRegC1)
        fprintf(' J1/J0=%-g \t Idata1/Idata0=%-f \t  IRegC1/IRegC0=%-f \n ',J1/J0,Idata1/Idata0,IRegC1/IRegC0)
        
        
        %input(' press return ')
        
        save JoptSave JoptVector
        
    end
    
    JoptVector=JoptVector(1:iteration+1);
    figure ; semilogy(0:numel(JoptVector)-1,JoptVector,'o-') ; title('J') ; xlabel('Iteration')
    
    
    
    % pause(30)
    
    % final calculation using last value for C
    %[u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen,C,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
    
    
    %figure(800) ; trisurf(TRIxy,x,y,Ctrue)    ;  title(' Ctrue')
    %%
    
    % 	%% Hessian calculations (seems very memory intensive, Hessian possible not very sparse)
    %
    % 	load TempSave S B h u v C m coordinates connectivity nip rho rhow CtrlVar lx ly dJdCCc K
    % 	[Bc]=calcBc(S,B,h,u,v,C,m,coordinates,connectivity,nip,rho,rhow,CtrlVar) ;
    % 	disp('Q')
    % 	Q=K\Bc;
    % 	disp('Hgn')
    % 	Hgn=Q'*Q; clear Q;
    % 	disp('dC')
    % 	dC=Hgn\dJdCCc;
    
    %%
    
    
    
    
    
    
    