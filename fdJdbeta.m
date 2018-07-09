function [dJdbeta,dJdb]=fdJdbeta(beta,beta_prior,b_prior,iModelType,iMisfitType,...
        BoundaryNodes,uModel,vModel,sModel,hModel,BModel,...
        uMeas,vMeas,wMeasInt,coordinates,connectivity,nip,...
        etaInt,gfint,AGlen,Luv,Luvrhs,n,m,alpha,rho,rhow,g,icount)
    %
    % calculates dJ/dbeta and dJ/db  where J=Imistfit + IReg using the adjoint method
    %
    
    global solutionphase
    persistent lambdau lambdaadjoint
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    Lambda=zeros(Nnodes,1) ; Mu=zeros(Nnodes,1) ;
    l=zeros(2*length(BoundaryNodes),1) ;% array for Lagrange multipliers
    
   C= ((sqrt(uModel.*uModel+vModel.*vModel)).^(1-m)).*beta.^(-2) ;
   
   % disp(' Step#1: Solve the forward model ')
    
   
   
    solutionphase=1 ; % solving state equation
    [uModel,vModel,lambdau,kv]=SSTREAM2d(sModel,hModel,uModel,vModel,coordinates,connectivity,nip,...
        etaInt,gfint,AGlen,C,Luv,Luvrhs,lambdau,n,m,alpha,rho,rhow,g,icount);
    
    %   [etaInt,xint,yint,exx,eyy,exy]=calcStrainRatesEtaInt(uModel,vModel,coordinates,connectivity,nip,AGlen,n);
    %   [wModel]=calcVerticalSurfaceVelocity(hModel,bModel,uModel,vModel,exx,eyy,xint,yint,coordinates,connectivity,nip);
      
    % solve Adjoint Equations
    
    solutionphase=2 ; % solving costate equation
      %[Lambda,Mu,l]=SolveAdjointEquationsIntegral(kv,Lambda,Mu,BoundaryNodes,l,uModel,vModel,uMeas,vMeas,coordinates,connectivity,nip,icount);
      [Lambda,Mu,lambdaadjoint]=SolveAdjointEquationsDiscrete(kv,Lambda,Mu,BoundaryNodes,l,coordinates,connectivity,nip,...
        hModel,BModel,uModel,vModel,uMeas,vMeas,wMeasInt,...
        Nnodes,icount,iMisfitType);
      
      
     [dImisfitdbeta]=betaGradientHilmar(uModel,vModel,Lambda,Mu,C,connectivity,coordinates,nip);
    
     [dImisfitdb]=bGradientHilmar(uModel,vModel,sModel,etaInt,rho,g,alpha,Lambda,Mu,connectivity,coordinates,nip);
     
     bModel=sModel-hModel; 
                                  
     [IReg,dIRegdbeta,dIRegdb]=ModelNorm(iModelType,beta,beta_prior,bModel,b_prior,coordinates,connectivity,nip);
                               
    
     x=coordinates(:,1); y=coordinates(:,2); DT = DelaunayTri(x,y); TRI=DT.Triangulation;
     figure(1001) ; trisurf(TRI,x,y,dImisfitdb) ;  title(' dImisfitb')
     figure(1002) ; trisurf(TRI,x,y,dImisfitdbeta) ;  title(' dImisfitdbeta')
     
    
    
    dJdbeta=dImisfitdbeta+dIRegdbeta;
    dJdb=dImisfitdb+dIRegdb;
    %dJdbeta=dImisfitdbeta; %+dIRegdbeta;
    
   %  x=coordinates(:,1); y=coordinates(:,2); DT = DelaunayTri(x,y); TRI=DT.Triangulation;
   %      figure(1001) ; trisurf(TRI,x,y,dJdbeta) ;  title(' dJdbeta')
    
   
   return
    
    % some stuff not needed by might come handy later
    
    
    
    
    %        Calculate beta gradient:
    %
    %        I've implemented three methods, Dough 1992, my version, and brute force calculation of directional
    %        derivatives of the matrix using complex analysis method
    %        All methods give same answers, but I continue to be confused as to why they differ
    %        so much from -2 beta_i (lambda_i u_i + mu_i v_i)
    %
    %         [betaGradDoug]=betaGradientDoug(uModel,vModel,Lambda,Mu,C,connectivity,coordinates,nip);
    %         figure(21) ; trisurf(TRI,x,y,betaGradDoug) ;  title(' betaGradDoug ')
    %         figure(22) ; trisurf(TRI,x,y,betaGradDoug./areas) ;  title(' betaGradDoug./areas ')
    
    %   The sensitivity of the cost function is suprising
    %          uModel=uModel*0+1;vModel=vModel*0+1;
    %          uMeas=uMeas*0+10*sin(2*pi*x/100e3);vMeas=vMeas*0+exp(-(y/100e3).^2);
    %          J0=CostFunction(sModel,uModel,vModel,wModel,bModel,BModel,...
    %              sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
    %              xMeas,yMeas,coordinates,connectivity,Nnodes,Nele,nip,nod);
    %
    %          dJ0db=zeros(Nnodes,1);
    %          for I=1:Nnodes
    %              uModelPert=uModel; vModelPert=vModel;
    %              uModelPert(I)=uModelPert(I)+1e-10;vModelPert(I)=vModelPert(I)+1e-10;
    %              dJ0db(I)=CostFunction(sModel,uModelPert,vModelPert,wModel,bModel,BModel,...
    %                  sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
    %                  xMeas,yMeas,coordinates,connectivity,Nnodes,Nele,nip,nod);
    %              dJ0db(I)=(dJ0db(I)-J0)/1e-10;
    %          end
    %
    %          save dJ dJ0db
    %          figure(101) ; trisurf(TRI,x,y,dJ0db) ;  title(' dJ0db')
    %
    
    %         Takes a lot of time, only good for verification purposes
    %         GradBeta=zeros(Nnodes,1);
    %         for I=1:Nnodes
    %             [DkvDbeta]=calcDkvDbeta(I,sModel,hModel,uModel,vModel,coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,C,m);
    %             GradBeta(I)=[Lambda;Mu]'*DkvDbeta*[uModel;vModel];
    %          end
    %
    %         figure(600) ; trisurf(TRI,x,y,GradBeta) ;  title(' GradBeta ')
    %          figure(601) ; trisurf(TRI,x,y,GradBeta-betaGradHilmar) ;  title(' GradBeta-betaGradHilmar ')
    %
    % simple method
    %         beta= sqrt(C.^(-1/m).* (sqrt(uModel.*uModel+vModel.*vModel)).^(1/m-1));
    %         sens= -2*beta.*(Lambda.*uModel+Mu.*vModel);
    %         figure(610) ; trisurf(TRI,x,y,sens) ;  title(' sens')
    %
    
    
end

