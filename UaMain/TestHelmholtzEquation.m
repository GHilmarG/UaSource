
%% Create mesh

MainFigure=10;    % Just the number of the main figure windwo
DoPlots=1;
UserVar=[];



% Create mesh
CtrlVar=Ua2D_DefaultParameters(); %
CtrlVar.MUA.MassMatrix=1;
CtrlVar.TriNodes=3;


xmin=-10e3; xmax=10e3 ; ymin=-10e3 ; ymax=10e3; ds=0.25e3; 
MeshBoundaryCoordinates=[xmin ymin ; xmin ymax ; xmax ymax ; xmax ymin];
CtrlVar.MeshSize=ds;
CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;

[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotMuaMesh(CtrlVar,MUA)

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

%%
close all

Test=9 ;
fPriorRHS=[];
switch Test
    case 1
        % test 1  : gives f=1
        a=1;
        b=1;
        c=1 ;
        d=0;
    case 2
        c=zeros(MUA.Nnodes,1);
        I=abs(x)<1 & abs(y)< 10 ;
        c(I)=1;
        d=0;
    case 3
        a=1 ;
        b=0.1 ;  % The smaller I make b, the smaller norm(f-c) becomes
        c=sin(2*pi*(x-xmin)/(xmax-xmin));
        d=0;
    case 4  % now f is equal to c irrespectivly of the value of scale 
        
        scale=1e4;
        c=scale*sin(2*pi*(x-xmin)/(xmax-xmin));
        b=0 ;
        a=scale*1;
        d=0;
    case 5  % now f is a smooth version of c
        
        a=1;
        b=1e6 ;  % b has the units l^{-2} compared to a unitless a
        c=DiracDelta(1/1e3,x,-1e3)+DiracDelta(1/1e3,x,1e3);
        d=0;
    case 6  % this will make the solution f equal to the fPrior, irrespectivly of the values of a and b
        
        fPrior=DiracDelta(1,x,-1)+DiracDelta(1,x,1);
        b=1+x*0;
        a=1+x*0;
        c=NaN+x*0; % does not affect the solution
        d=NaN;
        [UserVar,f,lambda,HEmatrix,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c,d);
        fPriorRHS=HEmatrix*fPrior;
        
    case 7  % inovation process spatial Gaussina white noise with unit variance
   

        fPrior=DiracDelta(1,x,-1)+DiracDelta(1,x,1);
        alpha=2 ; 
        rho=1e4; 
        dimention=2; 
        distance=linspace(100*eps,(xmax-xmin)/4,100) ;
        
        [CovMatern,nu,kappa,sigma2]=Matern(alpha,rho,dimention,distance);
        % kappa returned is a function of rho and alpha and dimention  
        figure ; plot(distance,CovMatern) ; title('Matern covariance') ; xlabel('distance')
        
        a=kappa^2+x*0;
        b=1+x*0;
        s0=1;
        c=s0*randn(MUA.Nnodes,1);
        d=0;
        [UserVar,f,lambda,HEmatrix,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c,d);

        
    case 8  % Example of fitting a curve
        
        % The idea is to mimic the minimisation process J= D + R
        % Where R uses a Matern covaricane and D is here just some other simple
        % contraint on f.
        %
        % R=(f-fPrior)'*HEmatrix*(f-fPrior)/2/err^2 ;
        % dRdf=HEmatrix*(f-fPrior)/err^2 ;
        % dRdf=0 -> HEmatrix*(f-fPrior) = 0 or f=fPrior
        %
        % cov(f-fPrior,f-fPrior)=Matern
        %
        % J = D  +  R
        %
        % for example D= || f ||
        % J =  f M f/2  +  (f-fPrior) HE (f-fPrior)/2
        % dJ =  M f + HE (f-fPrior)/err
        % dJ=0  (M+ME/err) f = HE fPrior/err
        %
        fPrior=DiracDelta(1/1000,x,-1000)+DiracDelta(1/1000,x,1000);
        
        alpha=2 ;
        rho=10e3;
        dimention=2;
        distance=linspace(100*eps,(xmax-xmin)/4,100) ;
        [CovMatern,nu,kappa,sigma2]=Matern(alpha,rho,dimention,distance);
        figure ; plot(distance,CovMatern) ; title('Matern covariance') ; xlabel('distance')
        
        
        a=kappa^2+x*0;
        b=1+x*0;
        c=NaN;
        d=NaN;
        
        [UserVar,f,lambda,HE,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c,d);
        errNoise=x*0+1;
        I=x>0 ; errNoise(I)=1e-6;
        Sigma=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1./errNoise.^2) ;
        f=(MUA.M+Sigma*HE)\ (Sigma*HE*fPrior);
        
        figure ; plot(x/1000,f,'.r') ; hold on ; plot(x/1000,fPrior,'.b') ;  plot(x/1000,x*0,'.g') ;
        legend('f','prior','data')
        
    case 9  % using correlated and uncorrelated noise covariance
        
        % Using two covariance matrices,  
        % K= Matern  + S I
        % R=(f-fPrior)'*HE*(f-rPrior) / 2 +   (f-fData)' S * M *(f-fData)/2 ;
        %
        % I must be M because I need to evaluate an intergral \int f I f = f M f
        %
        %
        % dRdf=HE * (f-fPrior)  + S M*(f-fData)  ; 
        % dRdf=0 -> 0 = HE * (f-fPrior) + S M*(f-fData)  ;
        %              (HE + S M ) f = HE fPrior + S M fData
        %                          f = (HE+ S M) \ ( HE fPrior + S M fData)
        %
        % Here the prior is set equal to the noise-free data.  And then some noise is
        % added to the data
        %
        %
        fData=DiracDelta(1/1000,x,-1000)+DiracDelta(1/1000,x,1000);
        fData=fData/max(fData);
        fData=x*0+1;
        
        % an example using discontinous prior
        I=x>0; fPrior=x*0; fPrior(I)=fData(I) ; 
        
        % an example using discontinous prior
        fPrior=2*atan(x/1e3)/pi ; % continuous prior
        
        errNoise=0.1 ; fData=fData+errNoise*randn(numel(x),1);
        errPrior=x*0+0.01 ;  I=x<-5e3 ; errPrior(I)=errPrior(I)*1000; 
        
        PrecisionNoise=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1./errNoise.^2) ; % this is the precision of the noise
        PrecisionPrior=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1./(errPrior.^2)) ; % this is the precision of the Prior
        
        alpha=2 ;
        rho=1e+4;
        
        dimention=2;
        distance=linspace(100*eps,(xmax-xmin)/4,100) ;
        [CovMatern,nu,kappa,sigma2]=Matern(alpha,rho,dimention,distance);
        figure ; plot(distance,CovMatern) ; title('Matern covariance') ; xlabel('distance')
        
        
        a=kappa^2+x*0;
        b=1+x*0;
        c=NaN;
        d=NaN;
        
        
        MScale=sigma2; 
     
        a=a*MScale; b=b*MScale; % have to figure out how to affect the variance
        [UserVar,f,lambda,HE,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c,d);
        errNoise=errNoise+x*0;
        %I=x>0  ; err(I)=1e2;
         
        
   
        
        M=MUA.M; 

        f = ( PrecisionPrior*HE+ PrecisionNoise* M) \ ( PrecisionPrior* HE * fPrior + PrecisionNoise* M * fData) ; 
        
        figure ; plot(x/1000,fData,'.g') ; hold on ; plot(x/1000,f,'or') ; 
        plot(x/1000,fPrior,'.b') ;  
        plot(x/1000,f-fPrior,'.k')
        legend('data','f','prior','f-fPrior')
        xlabel('x (km)')
        title(' f = (PP HE + PN M ) \\ (PP HE fPrior + PN N fData)')
        fprintf(' norm(f-fData)=%f \t \t  norm(f-fPrior)=%f \n ',norm(f-fData),norm(f-fPrior))
        fprintf(' all gott. \n')
end

%%
[UserVar,f,lambda]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c,d,fPriorRHS);

figure(1010); PlotMeshScalarVariable(CtrlVar,MUA,f) ; title('f')
figure(1015); PlotMeshScalarVariable(CtrlVar,MUA,c) ; title('c')

figure(1020) ; plot(x,f,'.') ; title('f')
figure(1025) ; plot(x,c,'.') ; title('c')

CtrlVar.PlotXYscale=1000;
figure ; Plot_sbB(CtrlVar,MUA,f,f,f,[],[],1e-3); title('f')


%%
Norm=(f-fPrior)'*HEmatrix*(f-fPrior)/2
% derivative is HEmatrix*(f-fPrior)
% the derivative is 0 when HEmatrix*(f-fPrior)=0 -> f=fPrior
% 




%% comparision between semi-variogram and matern covariance
dd=ds/2 ; nd=round((xmax-xmin)/4/(dd)) ;  % dd*nd < (xmax-xmin)/2
[sv,sv_dist,nav]=get_sv(x,y,f,dd,nd);
figure(1050) ; plot(sv_dist,sv,'or') ; title('semi-variogram') ; xlabel('dist (length)') ; ylabel('(length^2)') 


hold on
% normalize 
plot(distance,CovMatern,'b')   ; % this does not give a good fit as the sigma2 expression appears incorrect
CovMatern=CovMatern/CovMatern(1);
s2=max(sv);
plot(distance,s2*(1-CovMatern),'r')

%%
Ff=scatteredInterpolant(x,y,f);
[X,Y]=ndgrid(xmin:ds:xmax,ymin:ds:ymax) ;
fGrid=Ff(X,Y);
figure(1030) ; contourf(X,Y,fGrid) ; colorbar ; axis equal

%Corr=xcorr2(fGrid-mean(mean(fGrid))); figure(1040) ; contourf(Corr') ; axis equal
