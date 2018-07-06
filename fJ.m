
function J=fJ(q,beta_prior,beta_error,b_prior,b_error,M,iModelType,iMisfitType, ...
        sModel,uModel,vModel,BModel,sMeas,uMeas,vMeas,wMeasInt,bMeas,BMeas,...
        coordinates,connectivity,nip,etaInt,gfint,AGlen,Luv,Luvrhs,n,m,alpha,rho,rhow,g,icount,n1beta,n2beta,n1b,n2b)
   
    
    beta=q(n1beta:n2beta);
    bModel=q(n1b:n2b);
    hModel=sModel-bModel;
    % Calculates the function J=Imisfit+IReg
    % This requires solving the SSTREAM equation and Imisfit
    
    global solutionphase
    persistent lambdau
    
    %beta= sqrt(C.^(-1/m).* (sqrt(u.*u+v.*v)).^(1/m-1)) ;
  
    C= ((sqrt(uModel.*uModel+vModel.*vModel)).^(1-m)).*beta.^(-2) ; 
    
%      x=coordinates(:,1); y=coordinates(:,2); DT = DelaunayTri(x,y); TRI=DT.Triangulation;
%     figure(800) ; trisurf(TRI,x,y,hModel-(sMeas-bMeas)) ;  title(' hModel-(sMeas-bMeas)')
%     figure(801) ; trisurf(TRI,x,y,sModel-sMeas) ;  title(' sModel-sMeas')
%      
    
    solutionphase=1 ; % solving state equation
    [uModel,vModel,lambdau,kv]=SSTREAM2d(sModel,hModel,uModel,vModel,coordinates,connectivity,nip,...
        etaInt,gfint,AGlen,C,Luv,Luvrhs,lambdau,n,m,alpha,rho,rhow,g,icount);
    

   
    [wModelInt,B] = VertVelMatrixVector(hModel,bModel,uModel,vModel,coordinates,connectivity,nip);
    % an alternative would be
    %[etaInt,xint,yint,exx,eyy,exy]=calcStrainRatesEtaInt(uModel,vModel,coordinates,connectivity,nip,AGlen,n);
    %[w2,wint2]=calcVerticalSurfaceVelocity(hModel,bModel,uModel,vModel,exx,eyy,xint,yint,coordinates,connectivity,nip);
    % wModelInt=wint2(:)
    
    
    
    %-----]
    
    
    Imisfit=ModelDataMisfit(sModel,bModel,hModel,uModel,vModel,B,M,sMeas,uMeas,vMeas,wMeasInt,bMeas,BMeas,coordinates,connectivity,nip,iMisfitType);
    
    [IRegbeta,IRegb]=ModelNorm(iModelType,beta,beta_prior,beta_error,bModel,b_prior,b_error,coordinates,connectivity,nip);
    
    
    
    
    J=Imisfit+IRegbeta+IRegb;
    
    disp([' J : ',num2str(J),', Imisfit : ',num2str(Imisfit),', IRegbeta : ',num2str(IRegbeta),' IRegb : ',num2str(IRegb)])
    
    
end