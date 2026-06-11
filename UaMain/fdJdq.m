function [dJdq]=fdJdq(q,beta_prior,beta_error,b_prior,b_error,M,iModelType,iMisfitType,...
        BoundaryNodes,uModel,vModel,sModel,BModel,...
        uMeas,vMeas,wMeasInt,coordinates,connectivity,nip,...
        etaInt,gfint,AGlen,Luv,Luvrhs,n,m,alpha,rho,rhow,g,icount,n1beta,n2beta,n1b,n2b)
    %
    % calculates dJ/dbeta and dJ/db  where J=Imistfit + IReg using the adjoint method
    %
    
    PlotGradients=0;
    
    global solutionphase
    persistent lambdau lambdaadjoint
    
    beta=q(n1beta:n2beta);
    bModel=q(n1b:n2b);
    hModel=sModel-bModel;
    
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
        hModel,sModel,uModel,vModel,uMeas,vMeas,wMeasInt,...
        Nnodes,icount,iMisfitType,M);
    
    
    [dImisfitdbeta]=betaGradientHilmar(uModel,vModel,Lambda,Mu,C,connectivity,coordinates,nip);
    
   % [dImisfitdb]=bGradientHilmar(uModel,vModel,sModel,hModel,etaInt,gfint,rho,rhow,g,alpha,m,C,Lambda,Mu,connectivity,coordinates,nip);
    [dImisfitdb]=bGradientHilmarBC(uModel,vModel,sModel,hModel,etaInt,gfint,rho,rhow,g,alpha,m,C,Lambda,Mu,connectivity,coordinates,nip);
   
    
    
    
    
    
  
    
    
    bModel=sModel-hModel;
    [~,~,dIRegdbeta,dIRegdb]=ModelNorm(iModelType,beta,beta_prior,beta_error,bModel,b_prior,b_error,coordinates,connectivity,nip);
    
    
    
    if PlotGradients==1
        
        x=coordinates(:,1); y=coordinates(:,2); DT = DelaunayTri(x,y); TRI=DT.Triangulation;
        
        
        figure(1001) ; trisurf(TRI,x,y,dImisfitdb) ;  title(' dImisfit/db')
        figure(1010) ; trisurf(TRI,x,y,dImisfitdbeta) ;  title(' dImisfit/dbeta')
        
        figure(1020) ; trisurf(TRI,x,y,dIRegdbeta) ;  title(' dIReg/dbeta')
        figure(1030) ; trisurf(TRI,x,y,dIRegdb) ;  title(' dIReg/db')
        
    end
    
    % Gradient Factors
    dbetaFactor=0 ; dbFactor=1;
    dImisfitdbeta=dImisfitdbeta*dbetaFactor  ; dImisfitdb=dImisfitdb*dbFactor;
    
    
    dJdbeta=dImisfitdbeta+dIRegdbeta;
    dJdb=dImisfitdb+dIRegdb;
    
    dJdq=[dJdbeta;dJdb]; 
    
    disp([' norm(dIRegdbeta) : ',num2str(norm(dIRegdbeta)),' norm(dIRegdb) : ',num2str(norm(dIRegdb))])
    disp([' norm(dImisfitdbeta) : ',num2str(norm(dImisfitdbeta)),' norm(dImisfitdb) : ',num2str(norm(dImisfitdb))])
    disp([' norm(dJdbeta) : ',num2str(norm(dJdbeta)),' norm(dJdb) : ',num2str(norm(dJdb))])
    
    
end

