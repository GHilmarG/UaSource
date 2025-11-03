
function GroundingLineCurvePlots(FileName)
    
    
    %
    
    if nargin==0
        [FileName,PathName,FilterIndes]=uigetfile('*.mat');
    end
    
    if isequal(FileName,0) ; return ; end
    
    % plotting setting for exporting
    
    %%
   
    
    CtrlVar.GLtension=1e-11;
    
    
    figPosition=[50 50 650 650];
    lw=10; % line width for .png does not work, use pdf for export
    % set path

% locdir=pwd;
%    indsGHG=strfind(upper(locdir),'GHG');
%    addpath([locdir(1:indsGHG+2),'/my_mathlab_functions'],'-append')
    
    %
    %
    fprintf(' Loading %s ',FileName)
    load(FileName,'as','ab','u','v','coordinates','connectivity','nip','AGlen','n','CtrlVar','rho','rhow','g','h','s','S','B','b','DTxy','TRIxy','MeshBoundaryCoordinates','dhdt','C','m','n','time','CtrlVar')
    fprintf(' done \n ')
    
    colormap('HSV')
    CtrlVar.PlotXYscale=1000;
    if ~isfield(CtrlVar,'PlotXYscale')
        CtrlVar.PlotXYscale=1000;
    end
    
    
   [~,BoundaryEdgeCornerNodes]=FindBoundaryNodes(connectivity,coordinates);  xBoundary=coordinates(BoundaryEdgeCornerNodes,1);  yBoundary=coordinates(BoundaryEdgeCornerNodes,2);  
   
    %
    % can get rid of this once all restart files contain all these fields
    
    [DTxy,TRIxy,DTint,DTintTriInside,Xint,Yint,xint,yint,Iint]=TriangulationNodesIntegrationPoints(coordinates,connectivity,BoundaryEdgeCornerNodes,nip);
    
    %[DTxy,TRIxy,DTint,DTintTriInside,Xint,Yint,xint,yint,Iint]=TriangulationNodesIntegrationPoints(coordinates,connectivity,MeshBoundaryCoordinates,nip);
    
    % calculate strain rates and deviatoric stresses, get rid of this once restart files contain these fields
    
    [etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
    
    
    
    %save TestSave
    %error('sdfa')
    %
    %load TestSave
    
    %%
    
    
    GF = GL2d(B,S,h,rhow,rho,connectivity,CtrlVar);
    CtrlVar.GLresolutionWhenPlotting=1000;
    [xglc,yglc,nxGL,nyGL]=FindNiceLookingGLforPlottingPurposes('Contouring',DTxy,GF,CtrlVar,xBoundary,yBoundary);
     
    %[xglc,yglc,nxGL,nyGL]=FindNiceLookingGLforPlottingPurposes('Contouring',DTxy,GF,CtrlVar);
    
    figure ; plot(xglc/1000,yglc/1000) ; axis equal
    
    
    
    
    
    %% This must be done for each grounding line
    
    [qngl,qtgl,ThetaN,ThetaT,kappaN,kappaT,Ngl,Tgl,f,hgl,dgl,ugl,vgl,rhogl,txxgl,tyygl,txygl]=...
        GLquantities(CtrlVar,xglc,yglc,nxGL,nyGL,coordinates,DTxy,DTint,Iint,s,b,S,B,h,u,v,g,rho,rhow,txx,tyy,txy);
    
    
    %error('asdf')
    
    %%
    % Plot quantities of interest
    
    
    % guess a reasonable spatial range for mesh plots
    fxmin=10*floor(min(xglc/10)/CtrlVar.PlotXYscale); fxmax=ceil(max(xglc/10)/CtrlVar.PlotXYscale)*10;
    fymin=10*floor(min(yglc/10)/CtrlVar.PlotXYscale); fymax=ceil(max(yglc/10)/CtrlVar.PlotXYscale)*10;
    fymin=max([fymin min(coordinates(:,2))/CtrlVar.PlotXYscale]) ;
    fxmin=max([fxmin min(coordinates(:,1))/CtrlVar.PlotXYscale]) ;
    fymax=min([fymax max(coordinates(:,2))/CtrlVar.PlotXYscale]) ;
    fxmax=min([fxmax max(coordinates(:,1))/CtrlVar.PlotXYscale]) ;
    if (fxmax-fxmin)/(fymax-fymin)<4
        fxmin=fxmin-(fymax-fymin)/8 ; fxmax=fxmax+(fymax-fymin)/8 ;
    end
    %
    
    plot(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,'r')
    set(gcf,'Name','xGL yGL')   ;
    xlabel('x (km)') ; ylabel('y (km)')
    axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])
    
    
    figure ; color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,f*0+100,f,'LineWidth',lw) ;
    set(gcf,'Name','floating ratio (0 fully floating, 1 grounded)')   ;
    colorbar
    xlabel('x (km)') ; ylabel('y (km)')
    %%
    
    % q= (A (rho g)^(n+1) (1-rho/rhow)^n' / (4^n' C') ))^(1/(1+m)) theta^(n/(m'+1)) h^((m'+n+3)/(m+1))  ( Schoof notation)
    %
    % where tau_bx=C' |\mathbf{v}|^{m'-1} u and |tau| = C' |v|^m' (Schoof)
    % I use: u=C |tau|^{m-1} tau_bx    
    % and I have |v|= C |tau|^m
    % therfore:  |v|= C (C' |v|^m')^m = C C'^m  |v|^(m' m)
    %->  m'=1/m and C'=(1/C)^(1/m)
    %
    % % q= A (rho g)^(n+1) (1-rho/rhow)^n / (4^n C^(-1/m)) ) theta^(n/(1/m+1)) h^((1/m+n+3)/(1/m+1))  ( my notation)
    %
    % n=3 ; m=3;
    % AGlen=1e-14*1e9*365.25*24*60*60; C0=7.624e6^(-m)*1000^m*365.25*24*60*60;
    % rho=900 ; rhow=1000; g=9.81/1000;
    lw=10;
    A0=mean(AGlen) ; C0=mean(C);
    qSchoof=rhogl.* (A0*(rhogl*g).^(n+1).*(1-rhogl/rhow).^n/(4^n*C0^(-1/m))).^(1/(1/m+1)).*ThetaN.^(n/(1/m+1)).*hgl.^((1/m+n+3)/(1/m+1)) ; % my notation
    qSchoof=real(qSchoof);
    
    
    
    %%
    
    
    [pathstr,fname,ext]=fileparts(FileName);
    fname=[fname,'gcurve'];
    fprintf('saving %s \n',fname)
    save(fname,'xglc','yglc','qngl','qtgl','ThetaN','ThetaT','kappaN','kappaT','Ngl','Tgl','f','hgl','dgl','ugl','vgl','txxgl','tyygl','txygl','qSchoof')
    
    
    %%
    
    
    T=(qngl-qSchoof)./qngl;
    
    figure
    CtrlVar.PlotNodes=0; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
    figure ; set(gcf,'Position',figPosition)
    PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ;
    hold on;
    color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,qngl*0+100,T,'LineWidth',lw) ; colorbar
    set(gcf,'Name','(q-qSchoof)/q')   ;
    colorbar ; xlabel('x (km)') ; ylabel('y (km)')
    axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])
    
    %%
    
    fprintf(' Plotting ice flux normal and tangential to grounding line \n ')
    
    CtrlVar.PlotNodes=0; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
    figure ; set(gcf,'Position',figPosition)
    PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ;
    hold on; htag=color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,qngl*0+100,qngl/1e9,'LineWidth',lw) ;
    %title('Ice flux normal to grounding line (10^9 kg m^{-1} a^{-1})')   ;
    set(gcf,'Name','Ice flux normal to grounding line (10^9 kg m^{-1} a^{-1})')   ;
    colorbar ; xlabel('x (km)') ; ylabel('y (km)')
    axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])
    
    
    
    CtrlVar.PlotNodes=0; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
    figure ; set(gcf,'Position',figPosition)
    PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ;
    hold on; htag=color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,qtgl*0+100,qtgl/1e9,'LineWidth',lw) ;
    title('Ice flux tangential to grounding line (10^9 kg m^{-1} a^{-1})')   ;
    set(gcf,'Name','Ice flux tangential to grounding line (10^9 kg m^{-1} a^{-1})')   ;
    colorbar ; xlabel('x (km)') ; ylabel('y (km)')
    axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])
    
    
    %%
    
    fprintf(' Plotting normal and tangential buttressing ratio \n ')
    
    
    CtrlVar.PlotNodes=0; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
    figure ; set(gcf,'Position',figPosition)
    PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ; daspect([1 1 1]) ;
    hold on; htag=color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,kappaN*0+100,kappaN,'LineWidth',lw) ;
    
    set(gcf,'Name','\kappa_n')   ;
    colorbar ; xlabel('x (km)') ; ylabel('y (km)')
    axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])
    
    CtrlVar.PlotNodes=0; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
    figure ; set(gcf,'Position',figPosition)
    PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ; daspect([1 1 1]) ;
    hold on; color_line3(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,kappaT*0+100,kappaT,'LineWidth',lw) ;
    %title('\Theta_t')   ;
    set(gcf,'Name','\kappa_t')   ;
    colorbar ; xlabel('x (km)') ; ylabel('y (km)')
    axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])
    
    
    fprintf(' Plotting deviatoric stresses along grounding line \n ')
    CtrlVar.PlotNodes=1; CtrlVar.PlotNodesSymbol='.'; CtrlVar.PlotNodesSymbolSize=3;
    figure ; set(gcf,'Position',figPosition)
    %PlotFEmesh(coordinates/CtrlVar.PlotXYscale,connectivity,CtrlVar) ; daspect([1 1 1]) ;
    
    % grounding line profile
    %plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'g','LineWidth',lw)
    plot(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,'g','LineWidth',2)
    %plot(xglc/CtrlVar.PlotXYscale,yglc/CtrlVar.PlotXYscale,'g','LineWidth',2)
    
    hold on
    scale=0.05 ;
    
    ind=1:50:numel(xglc);
    
    % label
    xls=(fxmin+0.1*(fxmax-fxmin))*CtrlVar.PlotXYscale; yls=(fymax-0.1*(fymax-fymin))*CtrlVar.PlotXYscale;
    xlst=(fxmin+0.2*(fxmax-fxmin))*CtrlVar.PlotXYscale; ylst=(fymax-0.1*(fymax-fymin))*CtrlVar.PlotXYscale;
    tx=100 ; ty=-100;
    text(xlst/CtrlVar.PlotXYscale,ylst/CtrlVar.PlotXYscale,'(100 kPa)','HorizontalAlignment','center')
    
    xlu=(fxmin+0.07*(fxmax-fxmin))*CtrlVar.PlotXYscale; ylu=(fymax-0.2*(fymax-fymin))*CtrlVar.PlotXYscale;
    xlut=(fxmin+0.2*(fxmax-fxmin))*CtrlVar.PlotXYscale; ylut=(fymax-0.2*(fymax-fymin))*CtrlVar.PlotXYscale;
    us=500 ; vs=0;
    text(xlut/CtrlVar.PlotXYscale,ylut/CtrlVar.PlotXYscale,'(500 m/a)','HorizontalAlignment','center')
    
    
    plot_strains([xls ;xglc(ind)]/CtrlVar.PlotXYscale,[yls;yglc(ind)]/CtrlVar.PlotXYscale,[tx;txxgl(ind)],[0;txygl(ind)],[ty;tyygl(ind)],scale,1);
    quiver([xlu;xglc(ind)]/CtrlVar.PlotXYscale,[ylu;yglc(ind)]/CtrlVar.PlotXYscale,[us;ugl(ind)],[vs;vgl(ind)],'color','k')
    axis([fxmin fxmax fymin fymax]) ; daspect([1 1 1])
    
    
    
    %title('Deviatoric stresses') ;
    xlabel('x (km)') ; ylabel('y (km)')
    set(gcf,'Name','Deviatoric stresses')   ;
    slope0a=973.6686e3;  slope0b=1265.712e3;
    y=coordinates(:,2);
    hold on ; plot([slope0a/CtrlVar.PlotXYscale,slope0a/CtrlVar.PlotXYscale],[min(y) max(y)]/CtrlVar.PlotXYscale,'m')
    hold on ; plot([slope0b/CtrlVar.PlotXYscale,slope0b/CtrlVar.PlotXYscale],[min(y) max(y)]/CtrlVar.PlotXYscale,'m')
    
    %%
    
end

%%
