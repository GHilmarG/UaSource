function CompareResultsWithAnalyticalTransferFunctions(Experiment,coordinates,connectivity,u,v,s,b,S,B,time,dt,AGlen,C,n,m,rho,rhow,g,alpha,nip,DTxy,TRIxy,DTint,TRIint,CtrlVar)
    
    if ~isequal(n,1) ; return ; end
    fprintf(' Comparing with analytical solutions \n ')
    x=coordinates(:,1); y=coordinates(:,2); h=s-b;
    
    PlotTriMesh=1; PlotTriContour=0;
    
    
    eta=1/(2*AGlen(1));
    
    
    % calculate numerical vertical velocities
    %
    [etaInt,xint,yint,exx,eyy,exy,Eint,e]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
    %[wSurf,wSurfInt,wBedInt]=calcVerticalSurfaceVelocity(rho,rhow,h,S,B,b,u,v,exx,eyy,xint,yint,coordinates,connectivity,nip,CtrlVar);
    as=s*0 ; ab=s*0;
    [wSurf,wSurfInt]=calcVerticalSurfaceVelocity(rho,rhow,h,S,B,b,u,v,as,ab,exx,eyy,xint,yint,coordinates,connectivity,nip,CtrlVar);
    
    %[X,Y]=meshgrid(x,y);
    Nfft=1024;
    xgrid=linspace(min(x),max(x),Nfft); ygrid=linspace(min(y),max(y),Nfft);
    [X,Y]=meshgrid(xgrid,ygrid);
    nx=length(xgrid) ; ny=length(ygrid) ; dx=xgrid(2)-xgrid(1);dy=ygrid(2)-ygrid(1);
    
    
    
    %%
    [sGrid,bGrid,SGrid,BGrid,alpha]=DefineGeometry(Experiment,[X(:) Y(:)],CtrlVar);
    [CGrid,m]=DefineSlipperyDistribution(Experiment,[X(:) Y(:)]);
    
    
    CGrid=reshape(CGrid,size(X));
    sGrid=reshape(sGrid,size(X));
    bGrid=reshape(bGrid,size(X));
    BGrid=reshape(BGrid,size(X));
    SGrid=reshape(SGrid,size(X));
    C0=mean(CGrid(:));
    b0=mean(bGrid(:));
    s0=mean(sGrid(:));
    h0=s0-b0;
    
    ub=C0*(mean(rho)*g*h0*sin(alpha))^m;
    
    
    db=bGrid-b0;
    dc=CGrid-C0;
    ds=sGrid-s0;
    
    
    
    delta_s=ifft2(ds); delta_b=ifft2(db); delta_c=ifft2(dc);
    k=fftspace(nx,dx); k(1)=eps ;
    l=fftspace(ny,dy); l(1)=eps ;
    
    % multiply transfer functions with corresponding basal perturbations
    % do inverse fft and only keep the real part
    
    
    sAnalytical=real(fft2(SSTREAM_Tss_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_s)+...
        fft2(SSTREAM_Tsb_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_b)+...
        fft2(SSTREAM_Tsc_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_c));
    
    uAnalytical=real(fft2(SSTREAM_Tus_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_s)+...
        fft2(SSTREAM_Tub_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_b)+...
        fft2(SSTREAM_Tuc_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_c));
    
    vAnalytical=real(fft2(SSTREAM_Tvs_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_s)+...
        fft2(SSTREAM_Tvb_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_b)+...
        fft2(SSTREAM_Tvc_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_c));
    
    
    wAnalytical=real(fft2(SSTREAM_Tws_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_s)+...
        fft2(SSTREAM_Twb_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_b)+...
        fft2(SSTREAM_Twc_t_3d_m(k,l,time,alpha,h0,eta,C0,mean(rho),g,m).*delta_c));
    
    sAnalytical=sAnalytical+h0 ;
    uAnalytical=uAnalytical+ub ;
    
    
    
    % interpolate analytical values onto the numerical FE mesh
    sAna = interp2(X,Y,sAnalytical,y,x,'spline');
    uAna = interp2(X,Y,uAnalytical,y,x,'spline');
    vAna = interp2(X,Y,vAnalytical,y,x,'spline');
    wAna = interp2(X,Y,wAnalytical,y,x,'spline');
    
    
    
    if PlotTriContour==1
        
        left=min(x) ; right=max(x) ; down=min(y) ; up=max(y);
        
        
        figure ;
        subplot(3,1,1); tricontour(coordinates,TRIxy,wSurf,15) ;  title(' w (nodal values) ') ;xlabel('x') ; ylabel('y');
        axis([left right down up]) ; axis equal tight ; colorbar
        
        subplot(3,1,2); contour(X,Y,wAnalytical',15) ;  title(' wAnalytical ') ;xlabel('x') ; ylabel('y');
        axis([left right down up]) ; axis equal tight ; colorbar
        
        
        subplot(3,1,3); tricontour(coordinates,TRIxy,wSurf-wAna,9) ;  title(' w-wAnalytical ') ;xlabel('x') ; ylabel('y');
        axis([left right down up]) ; axis equal tight ; colorbar
        
        figure ;
        subplot(3,1,1); tricontour(coordinates,TRIxy,u,15) ;  title(' u ') ;xlabel('x') ; ylabel('y');
        axis([left right down up]) ; axis equal tight ; colorbar
        
        subplot(3,1,2); contour(X,Y,uAnalytical',15) ;  title(' uAnalytical ') ;xlabel('x') ; ylabel('y');
        axis([left right down up]) ; axis equal tight ; colorbar
        
        
        subplot(3,1,3); tricontour(coordinates,TRIxy,u-uAna,9) ;  title(' u-uAnalytical ') ;xlabel('x') ; ylabel('y');
        axis([left right down up]) ; axis equal tight ; colorbar
        
        
        
        figure ;
        subplot(3,1,1); tricontour(coordinates,TRIxy,v,15) ;  title(' v ') ;xlabel('x') ; ylabel('y');
        axis([left right down up]) ; axis equal tight ; colorbar
        
        subplot(3,1,2); contour(X,Y,vAnalytical',15) ;  title(' vAnalytical ') ;xlabel('x') ; ylabel('y');
        axis([left right down up]) ; axis equal tight ; colorbar
        
        
        subplot(3,1,3); tricontour(coordinates,TRIxy,v-vAna,9) ;  title(' v-vAnalytical ') ;xlabel('x') ; ylabel('y');
        axis([left right down up]) ; axis equal tight ; colorbar
        
        
    end
    
    
    wAnaInt = interp2(X,Y,wAnalytical,yint(:),xint(:),'spline');
    
    if PlotTriMesh==1
        
        
        figure;
        subplot(3,4,1) ; trimesh(TRIxy,x,y,s); title('s (numerical)') ; xlabel('x') ; ylabel('y')
        subplot(3,4,5); trimesh(TRIxy,x,y,sAna); title('s (analytical)') ; xlabel('x') ; ylabel('y')
        subplot(3,4,9); trimesh(TRIxy,x,y,s-sAna); title('s (numerical-analytical)') ; xlabel('x') ; ylabel('y')
        
        subplot(3,4,2); trimesh(TRIxy,x,y,u-mean(u)); title('du (numerical)') ; xlabel('x') ; ylabel('y')
        subplot(3,4,6); trimesh(TRIxy,x,y,uAna-mean(uAna)); title('du (analytical)') ; xlabel('x') ; ylabel('y')
        subplot(3,4,10); trimesh(TRIxy,x,y,u-uAna); title('u (numerical-analytical)') ; xlabel('x') ; ylabel('y')
        
        subplot(3,4,3); trimesh(TRIxy,x,y,v-mean(v)); title('dv (numerical)') ; xlabel('x') ; ylabel('y')
        subplot(3,4,7); trimesh(TRIxy,x,y,vAna-mean(vAna)); title('v (analytical)') ; xlabel('x') ; ylabel('y')
        subplot(3,4,11); trimesh(TRIxy,x,y,v-vAna); title('dv (numerical-analytical)') ; xlabel('x') ; ylabel('y')
        
        subplot(3,4,4); trimesh(TRIxy,x,y,wSurf-mean(wSurf)); title('dw (numerical)') ; xlabel('x') ; ylabel('y')
        subplot(3,4,8); trimesh(TRIxy,x,y,wAna-mean(wAna)); title('w (analytical)') ; xlabel('x') ; ylabel('y')
        subplot(3,4,12); trimesh(TRIxy,x,y,wSurf-wAna); title('dw (numerical-analytical)') ; xlabel('x') ; ylabel('y')
        
        
    end
    
    velnorm=max(uAna)-min(uAna);
    
    u=u-mean(u) ; uAna=uAna-mean(uAna);
    v=v-mean(v) ; vAna=vAna-mean(vAna);
    wSurf=wSurf-mean(wSurf) ; wAna=wAna-mean(wAna);
    
    diffs=norm(s-sAna)/sqrt(length(s))/(max(b)-min(b));
    diffu=norm(u-uAna)/sqrt(length(u))/velnorm;
    diffv=norm(v-vAna)/sqrt(length(v))/velnorm;
    diffwInt=norm(wSurfInt(:)-wAnaInt(:))/sqrt(length(wAnaInt))/velnorm;
    diffw=norm(wSurf-wAna)/sqrt(length(wAna))/velnorm;
    %             diffw=norm(wSurfInt(:)-wAna(:))/sqrt(length(wAna(:)));
    
    fprintf('        time: %-g \n ',time)
    disp(['        rms s: ',num2str(diffs)])
    disp(['         rms u: ',num2str(diffu)])
    disp(['         rms v: ',num2str(diffv)])
    disp(['   rms w (int): ',num2str(diffwInt)])
    disp(['  rms w (node): ',num2str(diffw)])
    %             disp([' rms w: ',num2str(diffw)])
    
    [Nele,nod]=size(connectivity);
    
    fid = fopen('CompareAnalyticalNumerical.txt', 'a+');
    fprintf(fid, '%i \t %i \t %g \t %g \t %g \t %g \t %i \n', Nele,nod,diffs,diffu,diffv,diffw,nip);
    fclose(fid);
    
    
    % view the contents
end

%   Gauss b pert     100x100km alpha=0.001; hmean=1000; ampl_b=100; sigma_bx=10000 ; sigma_by=10000;
%    Nele        nod       u rms     v rms    w rms
%    2048         6
%
%

