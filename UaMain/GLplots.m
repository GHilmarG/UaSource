
function GLplots(FileName)
%%       

if nargin==0
        [FileName,PathName,FilterIndes]=uigetfile('*.mat');
end

if isequal(FileName,0) ; return ; end


%% set path
locdir=pwd;

indsGHG=strfind(upper(locdir),'GHG');
addpath([locdir(1:indsGHG+2),'/my_mathlab_functions'],'-begin')

load(FileName,'as','ab','u','v','coordinates','connectivity','nip','AGlen','n','CtrlVar','rho','rhow','g','h','s','S','B','b','DTxy','TRIxy','MeshBoundaryCoordinates','dhdt','C','m','n','time','GFstart','CtrlVar')


CtrlVar.PlotDeviatoricStresses=0;
CtrlVar.calcDerivedGLquantitiesForGLeleOnly=0;
H0=1000;
CtrlVar.GLthreshold=0.5;
% calculate strain rates and deviatoric stresses

[etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
GF = GL2d(B,S,h,rhow,rho,connectivity,CtrlVar);
GLgeo=GLgeometry(connectivity,coordinates,GF,CtrlVar);
[Fx,Fy,Ngl,Tgl,Theta,Ffree,f,qgl,nx,ny]=...
    calcDerivedGLquantities(CtrlVar,s,b,h,S,B,u,v,rho,rhow,g,coordinates,connectivity,nip,GF,GLgeo,txx,tyy,txy);
    

if CtrlVar.PlotDeviatoricStresses==1
    CtrlVar.PlotNodes=0;
    figure ; PlotFEmesh(coordinates/H0,connectivity,CtrlVar) ; axis equal tight ;
    
    % grounding line profile
    plot(GLgeo(:,[3 4])'/H0,GLgeo(:,[5 6])'/H0,'g','LineWidth',2)
    
    hold on
    scale=0.1 ;
    
    xGLele=mean(xint,2);
    yGLele=mean(yint,2);
    txxgl=mean(txx,2);
    tyygl=mean(tyy,2);
    txygl=mean(txy,2);
    
        
    N=10;
    plot_strains(xGLele(GLgeo(1:N:end,1))/H0,yGLele(GLgeo(1:N:end,1))/H0,txxgl(GLgeo(1:N:end,1)),txygl(GLgeo(1:N:end,1)),tyygl(GLgeo(1:N:end,1)),scale,2);
    axis equal
    title('Deviatoric stresses')
    
    
end

figure
quiver(mean(xint,2)/H0,mean(yint,2)/H0,nx,ny)
hold on ; plot(GLgeo(:,[3 4])'/H0,GLgeo(:,[5 6])'/H0,'k','LineWidth',2) ; xlabel('(km)') ; ylabel('(km)') ; title('Normals to GL.node')
scale=1;  hold on ; quiver(GLgeo(:,7)/H0,GLgeo(:,8)/H0,GLgeo(:,9),GLgeo(:,10),scale,'color','r')


figure ;
subplot(4,1,1) ;PlotElementBasedQuantities(coordinates/H0,connectivity,Ffree) ; title('Ffree') ; colorbar ;  % caxis([0 250])
hold on ; plot(GLgeo(:,[3 4])'/H0,GLgeo(:,[5 6])'/H0,'k','LineWidth',2) ; xlabel('(km)') ; ylabel('(km)')

subplot(4,1,2) ; PlotElementBasedQuantities(coordinates/H0,connectivity,Ngl)   ; title('Ngl')   ; colorbar ; % caxis([0 250])
hold on ; plot(GLgeo(:,[3 4])'/H0,GLgeo(:,[5 6])'/H0,'k','LineWidth',2) ; xlabel('(km)') ; ylabel('(km)')

subplot(4,1,3) ; PlotElementBasedQuantities(coordinates/H0,connectivity,Tgl)   ; title('Tgl')   ; colorbar ; % caxis([0 250])
hold on ; plot(GLgeo(:,[3 4])'/H0,GLgeo(:,[5 6])'/H0,'k','LineWidth',2) ; xlabel('(km)') ; ylabel('(km)')

subplot(4,1,4) ; PlotElementBasedQuantities(coordinates/H0,connectivity,Theta)   ; title('\theta')   ; colorbar ;
hold on ; plot(GLgeo(:,[3 4])'/H0,GLgeo(:,[5 6])'/H0,'k','LineWidth',2) ; xlabel('(km)') ; ylabel('(km)')

figure ; PlotElementBasedQuantities(coordinates/H0,connectivity,Theta)   ; title('\theta')   ; colorbar ;
hold on ; plot(GLgeo(:,[3 4])'/H0,GLgeo(:,[5 6])'/H0,'k','LineWidth',2) ; xlabel('(km)') ; ylabel('(km)')


figure ; PlotElementBasedQuantities(coordinates/H0,connectivity,f)   ; title('floating ratio')   ; colorbar ;
hold on ; plot(GLgeo(:,[3 4])'/H0,GLgeo(:,[5 6])'/H0,'k','LineWidth',2) ; xlabel('(km)') ; ylabel('(km)')


figure ; plot(Theta(GLgeo(:,1)),'o') ; title(' \theta for GL elements') ; xlabel('(km)') ; ylabel('(km)')

figure ; plot3(GLgeo(:,7)/H0,GLgeo(:,8)/H0,Theta(GLgeo(:,1)),'o') ; title(' \theta for GL elements') ; xlabel('(km)') ; ylabel('(km)')


figure ;PlotElementBasedQuantities(coordinates/H0,connectivity(GLgeo(:,1),:),qgl/1e9)  ; title('Ice flux at grounding line (10^9 kg m^{-1} a^{-1})')   ; 
colorbar('southoutside') ; axis equal tight
hold on ; plot(GLgeo(:,[3 4])'/H0,GLgeo(:,[5 6])'/H0,'k','LineWidth',2) ; xlabel('(km)') ; ylabel('(km)')



%%
% % calculating a few values of interest for the Ex3a experiments only
% ind=abs(GLgeo(:,8)) < 10e3; 
% hEleGL=mean(hEle(GLgeo(ind,1)));
% ThetaGL=mean(Theta(GLgeo(ind,1)));
% 
% 
% fprintf('       GL min x pos \t \t \t GL flux \t  \t std   \t \t ice thick \t \t theta \n')
% fprintf('%20.10g \t %20.10g \t %20.10g \t %20.10g \t %20.10g \n',minGLx/1000,mean(qgl(ind))/1e9,std(qgl(ind))/1e9,hEleGL,ThetaGL)
% 




%%
