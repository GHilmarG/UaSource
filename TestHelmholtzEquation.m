
%% Create mesh

MainFigure=10;    % Just the number of the main figure windwo
DoPlots=1;
UserVar=[];



% Create mesh
CtrlVar=Ua2D_DefaultParameters(); %
CtrlVar.TriNodes=6;


xmin=-10e3; xmax=10e3 ; ymin=-10e3 ; ymax=10e3; ds=0.25e3; 
MeshBoundaryCoordinates=[xmin ymin ; xmin ymax ; xmax ymax ; xmax ymin];
CtrlVar.MeshSize=ds;
CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;

[UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates);
figure ;  PlotMuaMesh(CtrlVar,MUA)

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

%%
close all

Test=4 ;
fPriorRHS=[];
switch Test
    case 1
        % test 1  : gives f=1
        a=1;
        b=1;
        c=1 ;
    case 2
        c=zeros(MUA.Nnodes,1);
        I=abs(x)<1 & abs(y)< 10 ;
        c(I)=1;
    case 3
        a=1 ;
        b=0.1 ;  % The smaller I make b, the smaller norm(f-c) becomes
        c=sin(2*pi*(x-xmin)/(xmax-xmin));
        
    case 4  % now f is equal to c irrespectivly of the value of scale 
        
        scale=1e4;
        c=scale*sin(2*pi*(x-xmin)/(xmax-xmin));
        b=0 ;
        a=scale*1;
        
    case 5  % now f is a smooth version of c
        
        a=1;
        b=1e6 ;  % b has the unites l^{-2} compared to a unitless a
        c=DiracDelta(1/1e3,x,-1e3)+DiracDelta(1/1e3,x,1e3);
        
    case 6  % this will make the solution f equal to the fPrior, irrespectivly of the values of a and b
        
        fPrior=DiracDelta(1,x,-1)+DiracDelta(1,x,1);
        b=1+x*0;
        a=1+x*0;
        c=NaN+x*0; % does not affect the solution
        [UserVar,f,lambda,HEmatrix,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c);
        fPriorRHS=HEmatrix*fPrior;
        
    case 7  % inovation process spatial Gaussina white noise with unit variance
        
        fPrior=DiracDelta(1,x,-1)+DiracDelta(1,x,1);
        alpha=2 ; rho=1e3; dimention=2; distance=linspace(100*eps,(xmax-xmin)/4,100) ;
        [CovMatern,nu,kappa,sigma2]=Matern(sigma,alpha,rho,dimention,distance);
        %kappa=0.0028284;
        a=kappa^2+x*0;
        b=1+x*0;
        s0=1;
        c=s0*randn(MUA.Nnodes,1);
        
        [UserVar,f,lambda,HEmatrix,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c);
        
        
    case 8  % inovation process spatial Gaussina white noise with unit variance
        
        fPrior=DiracDelta(1,x,-1)+DiracDelta(1,x,1);
        
        alpha=2 ; rho=1e3; dimention=2; distance=linspace(100*eps,(xmax-xmin)/4,100) ;
        [CovMatern,nu,kappa,sigma2]=Matern(sigma,alpha,rho,dimention,distance);
        
        a=kappa^2+x*0;
        b=1+x*0;
        s0=1;
        c=s0*randn(MUA.Nnodes,1);
        
        [UserVar,f,lambda,HEmatrix,HErhs]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c);
end


[UserVar,f,lambda]=HelmholtzEquation(UserVar,CtrlVar,MUA,a,b,c,fPriorRHS);

figure(1010); PlotMeshScalarVariable(CtrlVar,MUA,f) ; title('f')
figure(1015); PlotMeshScalarVariable(CtrlVar,MUA,c) ; title('c')

figure(1020) ; plot(x,f,'.') ; title('f')
figure(1025) ; plot(x,c,'.') ; title('c')

CtrlVar.PlotXYscale=1000;
figure ; Plot_sbB(CtrlVar,MUA,f,f,f,[],[],1e-3); title('f')

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
