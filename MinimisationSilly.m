function [Cest,AGlenEst,u,v,JoptVector]=MinimisationSilly(sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
        s,S,B,h,u,v,coordinates,connectivity,Xint,Yint,xint,yint,Boundary,DTxy,TRIxy,DTint,TRIint,...
        nip,AGlen,C,Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
        n,m,alpha,rho,rhow,g,GF,CtrlVar,Itime,JoptVector);
    
    if CtrlVar.InfoLevelAdjoint>=10;
        fprintf(' minimisation method simple \n ')
    end
    
    
    
    dIdAdescent=AGlen*0; dIdA=AGlen*0;  % currently always conj. grad. is always reset on restart
    dIdCdescent=C*0; dIdC=C*0;  % c
    
    wMeasInt=Grid1toGrid2(DTxy,wMeas,Xint,Yint);
    x=coordinates(:,1); y=coordinates(:,2);
    
    iA=strfind(CtrlVar.AdjointGrad,'A'); iC=strfind(CtrlVar.AdjointGrad,'C');
    isAgrad=~isempty(iA);
    isCgrad=~isempty(iC); 
    
    
    upA=zeros(length(C),1)+CtrlVar.AGlenmax;  lowA=zeros(length(C),1)+CtrlVar.AGlenmin;
    upC=zeros(length(C),1)+CtrlVar.Cmax;  lowC=zeros(length(C),1)+CtrlVar.Cmin;
    
    if norm(C - kk_proj(C,upC,lowC)) > 0
        disp(' initial C iterate not feasible ');
        C=kk_proj(C,upC,lowC);
    end
    
    
    if norm(AGlen - kk_proj(AGlen,upA,lowA)) > 0
        disp(' initial AGlen iterate not feasible ');
        AGlen=kk_proj(AGlen,upA,lowA);
    end
    
    
    
    nIt=CtrlVar.MaxAdjointIterations;
    gamma_Test=1;
    
    
    if CtrlVar.AdjointRestart==0;
        JoptVector=zeros(nIt+1,4)+NaN; iJ=0;
    else
        iJ=size(JoptVector,1)-1;
        [Nx,Ny]=size(JoptVector) ; if Ny==3 ; JoptVector=[JoptVector , zeros(Nx,1)]; end 
        JoptVector=[JoptVector;zeros(nIt,4)+NaN];
    end
    
    for iteration=1:nIt
        ig=0; fVector=zeros(100,4)+NaN; gamma_Vector=zeros(100,1)+NaN;
        C0=C;
        AGlen0=AGlen;
        
        [u,v,lambdauv,K]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen0,C0,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
        [f0,Idata0,IRegC0,IRegAGlen0]=MisfitFunction(u,v,[],uMeas,vMeas,wMeasInt,C0,C_prior,AGlen0,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
        fMin=f0 ; IdataMin=Idata0 ; IRegCmin=IRegC0; IRegAGlenmin=IRegAGlen0;
        gamma_Min=0;
        iJ=iJ+1;  JoptVector(iJ,1)=f0; JoptVector(iJ,2)=Idata0; JoptVector(iJ,3)=IRegC0; JoptVector(iJ,4)=IRegAGlen0;
        
        ig=ig+1; fVector(ig,1)=f0; fVector(ig,2)=Idata0 ; fVector(ig,3)=IRegC0;  fVector(ig,4)=IRegAGlen0;   gamma_Vector(ig)= 0;
        [tbx,tby,tb] = CalcBasalTraction(u,v,C,m,GF,CtrlVar);
        
        if CtrlVar.InfoLevelAdjoint>=10000
            
            if CtrlVar.doplots==1;
                figure ; trisurf(TRIxy,x,y,tb,'EdgeColor','none')  ;  title(sprintf('Iteration %-i \t tb',iteration)) ; lightangle(-45,30) ; lighting phong ;
                xlabel('x') ; ylabel('y')
                figure ; trisurf(TRIxy,x,y,tbx,'EdgeColor','none') ;  title(sprintf('Iteration %-i \t tbx',iteration)) ; lightangle(-45,30) ; lighting phong ;
                xlabel('x') ; ylabel('y')
                figure ; trisurf(TRIxy,x,y,tby,'EdgeColor','none') ;  title(sprintf('Iteration %-i \t tby',iteration)) ; lightangle(-45,30) ; lighting phong ;
                
                speed=sqrt(u.*u+v.*v);
                figure; trisurf(TRIxy,x,y,speed,'EdgeColor','none')  ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t speed',iteration))
                xlabel('x') ; ylabel('y')
                figure(999) ; trisurf(TRIxy,x,y,C0,'EdgeColor','none')  ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t C0',iteration))
                xlabel('x') ; ylabel('y')
            end
        end
        
        %% a simple rough estimate of dIdata/dC
        %         RegParameter=mean(tb(tb>eps))/10;
        %         %Ctemp=(uMeas.*tbx./uError./uError+vMeas.*tby./vError./vError)./(tb.^(m-1).*(tbx.*tbx./uError./uError+tby.*tby./vError./vError)+RegParameter^m);
        %         Ctemp=(uMeas.*tbx+vMeas.*tby)./(tb.^(m-1).*(tbx.*tbx+tby.*tby)+RegParameter^m);
        %         dIdatadC=2*(Ctemp-C0);
        %         dIdatadC=EleAverageInterpolate(dIdatadC,coordinates,connectivity,[],CtrlVar);
        %        figure; trisurf(TRIxy,x,y,dIdatadC,'EdgeColor','none')  ;lightangle(-45,30) ; lighting phong ;
        %        title(sprintf('Iteration %-i \t dIdatadC (not normalized)',iteration))
        %%
        
        %% Adjoint method
        
        [dJdC,dJdAGlen,u,v,wint,lx,ly,dIdC,dIdAGlen]=AdjointGradientNR2d(s,S,B,h,u,v,uMeas,vMeas,wMeasInt,Cd,CC,CAGlen,C_prior,AGlen_prior,coordinates,connectivity,Boundary,nip,AGlen,C0,...
            Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,n,m,alpha,rho,rhow,g,GF,Itime,CtrlVar,iteration,DTxy,TRIxy,K);
        
        if CtrlVar.doplots==1;
            figure(10001); trisurf(TRIxy,x,y,lx,'EdgeColor','none')      ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t lx',iteration))
            xlabel('x') ; ylabel('y')
            figure(10002); trisurf(TRIxy,x,y,ly,'EdgeColor','none')      ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t ly',iteration))
            xlabel('x') ; ylabel('y')
        end
        
        
%        [etaInt,xint,yint,exx,eyy,exy,Eint,e]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
%        save TestAdjoint x y lx ly TRIxy u v h etaInt xint yint AGlen rho exx eyy exy Eint e
%        error('fdsa')
        
        if isCgrad ;   dIdatadC=-dIdC; else  dIdatadC=zeros(length(C),1) ; end
        if isAgrad ;   dIdatadA=-dIdAGlen; else  dIdatadA=zeros(length(AGlen),1) ; end
        
    
        
        
        %%
        
        
        %%  Testing if the direction dIdatadC is a desent direction, and correcting the magnitude of the gradient
        % I=I(C0+gamma dIdCest)  => |dIdgamma|=dIdC'*dIdCest
        % define kappa as dIdC=kappa dIdCest (assuming it is in same direction)
        % if I calculate dIdgamma numerically I find: dIdgamma=dIdC'*dIdCest=kappa*dIdCest'*dIdCest and that
        % kappa=dIdgamma/(dIdCest'*dIdCest)
        % When I calutate dIdgamma numerically I need a small perturbation, if the size of that step is gamma_eps then
        % the step size after scaling is gamma_eps/kappa.
        % After having found the min I can test if that step size is much smaller than the gamma that gives the minumum value
        % If not then I used too large step size in calculating the slope.
        %
        
        dIdC=zeros(length(C),1) ; dIdA=zeros(length(AGlen),1) ; 
        
        if isCgrad
            %% rescale C gradient
            gamma_eps=1e-3*norm(C0)/norm(dIdatadC); % small perturbation to C0
            gamma=gamma_eps;
            Ctest=kk_proj(C0+gamma*dIdatadC,upC,lowC);  % changing C in the direction of dIdatadC
            [u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlen0,Ctest,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
            [feps,Idataeps]=MisfitFunction(u,v,[],uMeas,vMeas,wMeasInt,Ctest,C_prior,AGlen0,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
            
            dIdgamma=(feps-f0)/gamma_eps;
            dIdatadgamma=(Idataeps-Idata0)/gamma_eps;
            if dIdatadgamma >=0
                warning('MinimisationSimple:desent','not a desent direction, i.e. going in the direction dIdatadC does not decrease the Idata cost function. ')
            end
            
            scale=abs(dIdatadgamma)/(dIdatadC'*dIdatadC); % this involves an approximation assuming that dIdCest is in same direction as dIdC
            dIdatadC=scale*dIdatadC;
            gamma_eps=gamma_eps/scale;


            dIdCreg=CC\(C0-C_prior);
            dIdCdescentlast=dIdCdescent;
            dIdCdescent=-dIdatadC+CtrlVar.isRegC*dIdCreg;
            
                     
            if CtrlVar.AdjointConjugatedGradients==1 && iteration>1
                dIdClast=dIdC;
                [dIdC,ConjGradAngle,teta]=NewConjugatedGrad(dIdCdescent,dIdCdescentlast,dIdClast);
                fprintf(' ConjGradAngleC=%-g degrees, \t teta=%-g \n',ConjGradAngle,teta)
            else
                dIdC=-dIdCdescent; %  dIdatadA-CtrlVar.isRegAGlen*dIdAreg;
            end
            
            if CtrlVar.doplots==1;
            figure(1); trisurf(TRIxy,x,y,C0-C_prior,'EdgeColor','none')  ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t C0-C_prior',iteration))
            xlabel('x') ; ylabel('y')
            figure(2); trisurf(TRIxy,x,y,dIdatadC,'EdgeColor','none')  ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t dIdatadC',iteration))
            xlabel('x') ; ylabel('y')
            figure(3); trisurf(TRIxy,x,y,-dIdCreg,'EdgeColor','none')  ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t dIdCreg',iteration))
            xlabel('x') ; ylabel('y')
            figure(4); trisurf(TRIxy,x,y,dIdC,'EdgeColor','none')      ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t dIdatadC+dIdCreg',iteration))
            xlabel('x') ; ylabel('y')
            end
            
        end
        
         if isAgrad
            %% rescale AGlen gradient
            gamma_eps=1e-4*norm(AGlen0)/norm(dIdatadA); % small perturbation to AGlen0
            gamma=gamma_eps;
            AGlentest=kk_proj(AGlen0+gamma*dIdatadA,upA,lowA);  % changing AGlen in the direction of dIdatadA
            [u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlentest,C0,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
            [feps,Idataeps]=MisfitFunction(u,v,[],uMeas,vMeas,wMeasInt,C0,C_prior,AGlentest,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
            
            dIdatadgamma=(Idataeps-Idata0)/gamma_eps;
            if dIdatadgamma >=0
                warning('MinimisationSimple:desent','not a desent direction, i.e. going in the direction dIdatadA does not decrease the Idata cost function. ')
            end
            
            scale=abs(dIdatadgamma)/(dIdatadA'*dIdatadA); % this involves an approximation 
            fprintf(' scale for A is %-g \n ',scale)
            dIdatadA=scale*dIdatadA;
            gamma_eps=gamma_eps/scale;


            dIdAreg=CAGlen\(AGlen0-AGlen_prior);
            
            dIdAdescentlast=dIdAdescent;
            dIdAdescent=-dIdatadA+CtrlVar.isRegAGlen*dIdAreg;
                
            if CtrlVar.AdjointConjugatedGradients==1 && iteration>1
                dIdAlast=dIdA;
                [dIdA,ConjGradAngle,teta]=NewConjugatedGrad(dIdAdescent,dIdAdescentlast,dIdAlast);
                fprintf(' ConjGradAngleA=%-g degrees, \t teta=%-g \n',ConjGradAngle,teta)
            else
                dIdA=-dIdAdescent; %  dIdatadA-CtrlVar.isRegAGlen*dIdAreg;
            end
            
            if CtrlVar.doplots==1;
                figure(10); trisurf(TRIxy,x,y,AGlen0-AGlen_prior,'EdgeColor','none')  ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t AGlen0-AGlen_prior',iteration))
                xlabel('x') ; ylabel('y')
                figure(11); trisurf(TRIxy,x,y,dIdatadA,'EdgeColor','none')  ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t dIdatadA',iteration))
                xlabel('x') ; ylabel('y')
                figure(12); trisurf(TRIxy,x,y,-dIdAreg,'EdgeColor','none')  ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t dIdAreg',iteration))
                xlabel('x') ; ylabel('y')
                figure(13); trisurf(TRIxy,x,y,dIdA,'EdgeColor','none')      ;lightangle(-45,30) ; lighting phong ;  title(sprintf('Iteration %-i \t dIdatadA+dIdAreg',iteration))
                xlabel('x') ; ylabel('y')
            end
            
        end

        
        fprintf('norm(C0-C_prior)=%-g , \t',norm(C0-C_prior)) ; fprintf('norm(AGlen0-AGlen_prior)=%-g \n ',norm(AGlen0-AGlen_prior))
        
        
        
        %% simple line search
        
         Slope0=((dIdC+dIdA)'*(dIdC+dIdA));
         gamma_MinEstimate=0.5*f0/Slope0;
         % close to zero f=f0-Slope0*gamma,  aiming at reduction by kappa: kappa f0=f0-Slope gamma => gamma= f0 (1-kappa)/Slope
        
        
        gamma_a=0 ; gamma_c=gamma_MinEstimate; gamma_b=gamma_c/2 ; fa=f0 ; gamma_Eps=gamma_MinEstimate/1000;
        
        gamma=gamma_c;
        
        Ctest=kk_proj(C0+gamma*dIdC,upC,lowC);
        AGlentest=kk_proj(AGlen0+gamma*dIdA,upA,lowA);
        [u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlentest,Ctest,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
        [fc,Idata,IRegC,IRegAGlen]=MisfitFunction(u,v,[],uMeas,vMeas,wMeasInt,Ctest,C_prior,AGlentest,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
        
        ig=ig+1; fVector(ig,1)=fc; fVector(ig,2)=Idata ; fVector(ig,3)=IRegC;  fVector(ig,4)=IRegAGlen; gamma_Vector(ig)= gamma;
        iarm=0 ; iarmmax=15; iarmmin=0;
        target=0.75*f0;
        
        
        
        if  fc<=target
            fMin=fc ; IdataMin=Idata ; IRegCmin=IRegC;  IRegAGlenmin=IRegAGlen; gamma_Min=gamma;
            if CtrlVar.InfoLevelAdjoint>=10
                fprintf(' Initial step accepted. Misfit reduced from %-g to  %-g \t ratio fMin/f0=%-g \n ',f0,fMin,fMin/f0)
            end
        else
            if CtrlVar.InfoLevelAdjoint>=10
                fprintf(' Initial step not accepted entering line search. f0=%-g \t fc=%-g \t ratio fc/f0=%-g \n ',f0,fc,fc/f0)
            end
            
            iFminTry=1; iFminTryMax=5;
            while (fMin >  target && iarm<=iarmmax && iFminTry <= iFminTryMax) || iarm<iarmmin 
                iarm=iarm+1; iFminTry=iFminTry+1;
                gamma=gamma_b;
                Ctest=kk_proj(C0+gamma*dIdC,upC,lowC);
                AGlentest=kk_proj(AGlen0+gamma*dIdA,upA,lowA);
                [u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlentest,Ctest,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
                [fb,Idata,IRegC,IRegAGlen]=MisfitFunction(u,v,[],uMeas,vMeas,wMeasInt,Ctest,C_prior,AGlentest,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
                ig=ig+1; fVector(ig,1)=fb; fVector(ig,2)=Idata ; fVector(ig,3)=IRegC;  fVector(ig,4)=IRegAGlen; gamma_Vector(ig)= gamma;
                
                if fb < fMin ; fMin=fb ;  IdataMin=Idata ; IRegCmin=IRegC; IRegAGlenmin=IRegAGlen;  gamma_Min=gamma_b ; iFminTry=0 ; end
                
                if CtrlVar.InfoLevelAdjoint>=10;
                    fprintf(' It %-i \t iarm %-i \t f0=%-g \t fa=%-g \t fb=%-g \t fc=%-g \t fmin=%-g \t fb/f0=%-g \t fmin/f0=%-g \n',iteration,iarm,f0,fa,fb,fc,fMin,fb/f0,fMin/f0)
                    fprintf(' \t  \t \t  a=%-g \t  b=%-g \t  c=%-g  \t gamma_Min=%-g \t gamma_MinEstimate=%-g \n',gamma_a,gamma_b,gamma_c,gamma_Min,gamma_MinEstimate)
                end
                
                gamma_TestOld=gamma_Test;
                [gamma_Test,ParStatus ] = parabolamin(gamma_a,gamma_b,gamma_c,fa,fb,fc);
                fprintf(' Iteration %-i \t iarm=%-i : \t fitting a parabola gives gamma=%-g \n ',iteration,iarm,gamma_Test)
                
                if ParStatus==1  % parabolic fit did not give an acceptable minimum
                    if fc < fb && fb < fa      % decreasing values
                        gamma_Test=5*gamma_c;    % extrapolate
                    elseif fc > fb && fb > fa  % increasing values
                        gamma_Test=(gamma_a+gamma_b)/5;
                    else                       %
                        gamma_Test=(gamma_a+gamma_b)/3;
                    end
                end
                
                [gamma_Dist,imin]=min(abs(gamma_Test-gamma_Vector(2:ig-1)));
                if gamma_Dist < gamma_Eps
                    if fMin>0.999*f0
                        gamma_Eps=gamma_Eps/100;
                    else
                        fprintf(' breaking out of linesearch because new value (%-g) so close to a prevous one (%-g) \n',gamma_Test,gamma_Vector(1+imin))
                        iarm=iarm-1;
                        break
                    end
                end
                
                if gamma_Test<0 ; gamma_Test=(gamma_a+gamma_b)/10; end
                %if gamma_Min > 0.8*gamma_b ; gamma_Min=0.8*gamma_b ; elseif gamma_Min < 0.1*gamma_b ; gamma_Min=0.1*gamma_b; end
                
                if gamma_Test < gamma_b
                    gamma_c=gamma_b ; fc=fb; gamma_b=gamma_Test ; fa=f0 ; gamma_a=0;
                elseif gamma_Test >= gamma_b && gamma_Test< gamma_c
                    gamma_a=gamma_b ; fa=fb ; gamma_b=gamma_Test ;
                else  % extrapolation step
                    if gamma_Test> 5*gamma_c ; gamma_Test=5*gamma_c ; end
                    gamma_a=gamma_b ; fa=fb ; gamma_b=gamma_Test ; gamma_c=2*gamma_Test;
                    
                    gamma=gamma_c;
                    Ctest=kk_proj(C0+gamma*dIdC,upC,lowC);
                    AGlentest=kk_proj(AGlen0+gamma*dIdA,upA,lowA);
                    [u,v]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlentest,Ctest,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
                    [fc,Idata,IRegC,IRegAGlen]=MisfitFunction(u,v,[],uMeas,vMeas,wMeasInt,Ctest,C_prior,AGlentest,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
                    ig=ig+1; fVector(ig,1)=fc; fVector(ig,2)=Idata ; fVector(ig,3)=IRegC;  fVector(ig,4)=IRegAGlen; gamma_Vector(ig)= gamma;
                end
            end
            if ~(iarm<=iarmmax) ; fprintf(' Line search exit due to maximum number of iterations reached. \n') ; end
            if ~(iFminTry<=iFminTryMax) ; fprintf(' Exiting line search because no improvement found for the last %-i iterations. \n',iFminTry) ; end

           
        end
        
        if CtrlVar.InfoLevelAdjoint>=100 && CtrlVar.doplots==1;
            %% plotting
            nLS=5;
            
            gamma_StepVector=gamma_Min*[-1/3,1/3,2/3,4/3,5/3];
            fPlot=zeros(nLS,1) ; Idata=zeros(nLS,1) ; IRegC=zeros(nLS,1); IRegAGlen=zeros(nLS,1);
            parfor JJ=1:nLS
                
                gamma=gamma_StepVector(JJ);
                Ctest=kk_proj(C0+gamma*dIdC,upC,lowC);
                AGlentest=kk_proj(AGlen0+gamma*dIdA,upA,lowA);
                [uTemp,vTemp]=SSTREAM2dNR(s,S,B,h,u,v,coordinates,connectivity,Boundary,nip,AGlentest,Ctest,Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
                [fPlot(JJ),Idata(JJ),IRegC(JJ),IRegAGlen(JJ)]=MisfitFunction(uTemp,vTemp,[],uMeas,vMeas,wMeasInt,Ctest,C_prior,AGlentest,AGlen_prior,Cd,CAGlen,CC,coordinates,connectivity,nip,GF,CtrlVar);
                
            end
            
            fVector(ig+1:ig+nLS,1)=fPlot; fVector(ig+1:ig+nLS,2)=Idata ; fVector(ig+1:ig+nLS,3)=IRegC;  fVector(ig+1:ig+nLS,4)=IRegAGlen; gamma_Vector(ig+1:ig+nLS)= gamma_StepVector;
            
            ind=find(~isnan(gamma_Vector)) ; gamma_Vector=gamma_Vector(ind) ;
            
            temp=fVector;
            fVector=zeros(length(ind),4);
            fVector(:,1)=temp(ind,1); fVector(:,2)=temp(ind,2); fVector(:,3)=temp(ind,3); fVector(:,4)=temp(ind,4);
            [gamma_Vector,ind]=sort(gamma_Vector) ;
            fVector(:,1)=fVector(ind,1); fVector(:,2)=fVector(ind,2);fVector(:,3)=fVector(ind,3);fVector(:,4)=fVector(ind,4);
            figure(777); hold off ;
            plot(gamma_Vector,fVector(:,1),'-+b')
            hold on
            plot(gamma_Vector,fVector(:,2),'-og')
            plot(gamma_Vector,fVector(:,3),'-xr')
            plot(gamma_Vector,fVector(:,4),'-xc')
            legend('J','DataMisfit','SystemNormC','SystemNormA')
            hold on ; plot(gamma_Min,fMin,'*r') ; title(sprintf(' Iteration %-i ',iteration)) ; xlabel('\gamma') ; ylabel(' J ')
            plot([0,gamma_eps],[f0,f0-Slope0*gamma_eps],'r')
        end
        %%
        fprintf(' At the end of iteration %-i \t J=%-g, \t Jdata=%-g, \t IRegC=%-g, \t IRegAGlen=%-g, with fmin/f0=%-g \n',iteration,fMin,IdataMin,IRegCmin,IRegAGlenmin,fMin/f0)
        
        C=kk_proj(C0+gamma_Min*dIdC,upC,lowC);
        AGlen=kk_proj(AGlen0+gamma_Min*dIdA,upA,lowA);
        Cest=C ; AGlenEst=AGlen;
        
        fprintf(' \t  \t \t  a=%-g \t  b=%-g \t  c=%-g  \t gamma_Min=%-g \t gamma_MinEstimate=%-g \n',gamma_a,gamma_b,gamma_c,gamma_Min,gamma_MinEstimate)

    end
    iJ=iJ+1;  JoptVector(iJ,1)=fMin; JoptVector(iJ,2)=IdataMin; JoptVector(iJ,3)=IRegCmin;
    
    if CtrlVar.doplots==1;
        figure(100) ; trisurf(TRIxy,x,y,tb,'EdgeColor','none')  ;  title('tb')    ; xlabel('x') ; ylabel('y')
        figure(101) ; trisurf(TRIxy,x,y,tbx,'EdgeColor','none') ;  title('tbx')  ; xlabel('x') ; ylabel('y')
        figure(102) ; trisurf(TRIxy,x,y,tby,'EdgeColor','none') ;  title('tby')  ; xlabel('x') ; ylabel('y')
        
        
        figure(103) ; trisurf(TRIxy,x,y,dIdC.*GF.node,'EdgeColor','none')  ;  title('dIdC') ;lightangle(-45,30) ; lighting phong ;  xlabel('x') ; ylabel('y')
        figure(104) ; trisurf(TRIxy,x,y,dIdA,'EdgeColor','none')  ;  title('dIdA') ;lightangle(-45,30) ; lighting phong ;  xlabel('x') ; ylabel('y')
        
        figure(105) ; trisurf(TRIxy,x,y,u,'EdgeColor','none') ;  title('u') ; lightangle(-45,30) ; lighting phong ;  xlabel('x') ; ylabel('y')
        figure(106) ; trisurf(TRIxy,x,y,v,'EdgeColor','none') ;  title('v') ; lightangle(-45,30) ; lighting phong ;  xlabel('x') ; ylabel('y')
        
        figure(107) ; trisurf(TRIxy,x,y,u-uMeas,'EdgeColor','none') ;  title('u-uMeas') ; lightangle(-45,30) ; lighting phong ;  colorbar  ; xlabel('x') ; ylabel('y')
        figure(108) ; trisurf(TRIxy,x,y,v-vMeas,'EdgeColor','none') ;  title('v-vMeas') ; lightangle(-45,30) ; lighting phong ;  colorbar  ; xlabel('x') ; ylabel('y')
        
        figure(120) ; trisurf(TRIxy,x,y,Cest,'EdgeColor','none') ;  title('Cest') ; lightangle(-45,30) ; lighting phong ; colorbar  ; xlabel('x') ; ylabel('y')
        figure(121) ; trisurf(TRIxy,x,y,AGlenEst,'EdgeColor','none') ;  title('AGlenest') ; lightangle(-45,30) ; lighting phong ; colorbar  ; xlabel('x') ; ylabel('y')
        
        figure ;
        semilogy(0:iJ-1,JoptVector(1:iJ,1),'-+b') ; hold on
        semilogy(0:iJ-1,JoptVector(1:iJ,2),'-og') ; hold on
        semilogy(0:iJ-1,JoptVector(1:iJ,3),'-xr') ; hold on
        legend('J','DataMisfit','SystemNorm')
    end
    
    
    save SimpleResult Cest coordinates connectivity TRIxy dIdC u v
    
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
