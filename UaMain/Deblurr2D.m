

function [db,dc] = Deblurr2D(time,...
        sModel,uModel,vModel,bModel,BModel,...
        sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,xMeas,yMeas,...
        Experiment,...
        coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,AGlen,C,...
        Luv,Luvrhs,lambdau,n,m,alpha,rho,rhow,g,Itime)
    
   bModelStart=bModel ;  CModelStart=C ; 
    
%     load Deblurr2D time ...
%         sModel uModel vModel bModel BModel ...
%         sMeas uMeas vMeas wMeas bMeas BMeas xMeas yMeas ...
%         Experiment ...
%         coordinates connectivity Nnodes Nele nip nod etaInt gfint AGlen C ...
%         Luv Luvrhs lambdau n m alpha rho rhow g Itime bModelStart CModelStart
    
    
    disp('Deblur2D')
    
    for It=1:1
        
        
        hModel=sModel-bModel;
        
        [uModel,vModel,exx,eyy,exy,etaInt,beta2,xint,yint,kv,lambdau]=...
            SSS2diagnosticBCsparse(sModel,hModel,uModel,vModel,coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,AGlen,C,...
            Luv,Luvrhs,lambdau,n,m,alpha,rho,rhow,g,Itime);
        
        [wModel]=calcVerticalSurfaceVelocity(hModel,bModel,uModel,vModel,exx,eyy,xint,yint,coordinates,connectivity,nip);
        
      
        
        seps=1e10; ueps=1e-5 ; veps=ueps; weps=ueps;
        gb=1e30; gc=1e30;
        gs=0; gu=1/ueps^2; gv=1/veps^2; gw=1/weps^2 ;
        
        % for testing purposes synthetic forward data can be generated
        GenForwardData=0;
        
        
        
        xModel=coordinates(:,1) ; yModel=coordinates(:,2);
        Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
        
        
        % 1) rotate by the angle theta in anticlockwise direction round z
        %    after that the x' axis is along directly of steepest slope (dip) and y' gives the strike direction
        %    Rotation Matrix= [ cos(theta)  sin(theta)  0 ]
        %                     [-sin(theta)  cos(theta)  0 ]
        %                     [   0             0       1 ]
        % 2) then rotate the angle alpha0 around y' axis
        %    Rotation Matrix= [ cos(alpha0)  0   -sin(alpha0) ]
        %                     [    0        1      0        ]
        %                     [  sin(alpha0) 0    cos(alpha0) ]
        %
        % The problem with the second rotation is that as I rotate (x_i,y_i,s_i) I get different x,y values than
        % when I rotate (x_i,y_i,b_i), and values are not longer on the same grid.
        % I just don't do this as the rotation is very small. All I do is detrending.
        %
        % fit best plane to the surface and estimate: H0,alpha0, C0, eta
        
        %s-a0-ax*x-ay*y=0;
        sol=-[ones(Nnodes,1) xMeas yMeas] \ sMeas  ;
        ax=sol(2) ; ay=sol(3);
        theta=atan2(ay,ax); alpha0=atan(sqrt(ax^2+ay^2));
        disp([' theta : ',num2str(theta),' alpha0 : ',num2str(alpha0)])
        
        sMeasMean=mean(sMeas);
        H=mean(sMeas-bMeas);
        bMeasMean=sMeasMean-H;
        
        rho=917; g=9.81/1000;
        U=mean(uMeas*cos(theta)+vMeas*sin(theta));
        eta=mean(etaInt(:));
        C0=mean(C);
        
        % rotate velocity by theta around z axis
        
        uData=uMeas*cos(theta)+vMeas*sin(theta);
        vData=-uMeas*sin(theta)+vMeas*cos(theta);
        wData=wMeas;
        
        uCalc=uModel*cos(theta)+vModel*sin(theta);
        vCalc=-uModel*sin(theta)+vModel*cos(theta);
        wCalc=wModel;
        
        % rotate velocity by -alpha0 around y axis
        
        uMeas=uData ;  vMeas=vData ; wMeas=wData ;
        uModel=uCalc;  vModel=vCalc; wModel=wCalc;
        
        uData=uMeas*cos(alpha0)-wMeas*sin(alpha0);
        wData=uMeas*sin(alpha0)+wMeas*cos(alpha0);
        vData=vMeas;
        
        uCalc=uModel*cos(alpha0)-wModel*sin(alpha0);
        wCalc=uModel*sin(alpha0)+wModel*cos(alpha0);
        vCalc=vModel;
        
        % remove mean slope from s, b and B
        
        sData=sMeas+ax*xMeas+ay*yMeas;      bData=bMeas+ax*xMeas+ay*yMeas;      BData=BMeas+ax*xMeas+ay*yMeas;
        sCalc=sModel+ax*xModel+ay*yModel;   bCalc=bModel+ax*xModel+ay*yModel;   BCalc=BModel+ax*xModel+ay*yModel;
        
        
        %DT = DelaunayTri(xMeas,yMeas); TRI=DT.Triangulation;
        %figure ; trisurf(TRI,xMeas,yMeas,sData) ;  title('sData') ;xlabel('x') ; ylabel('y')
        
        
        % now create a 2D grid and interpolate all variables on the grid,
        % only define values in the convex set, elsewhere but in mean values
        
        xrotModel=xModel*cos(theta)+yModel*sin(theta); yrotModel=-xModel*sin(theta)+yModel*cos(theta);
        xrotMeas=xMeas*cos(theta)+yMeas*sin(theta); yrotMeas=-xMeas*sin(theta)+yMeas*cos(theta);
        
        xl=max(xrotModel)-min(xrotModel); yl=max(yrotModel)-min(yrotModel);
        xmin=min(xrotModel)-xl/2 ; xmax=max(xrotModel)+xl/2 ;
        ymin=min(yrotModel)-yl/2 ; ymax=max(yrotModel)+yl/2 ;
        
        dx=H ; dy=H ;
        nx=2^nextpow2(xl/dx);
        ny=2^nextpow2(yl/dy) ;
        
        xfft=linspace(xmin,xmax,nx); yfft=linspace(ymin,ymax,ny); dx=xfft(2)-xfft(1); dy=yfft(2)-yfft(1);
        [X,Y]=ndgrid(xfft,yfft);
        
        % interpolate measured quantities on a grid
        DT = DelaunayTri(xrotModel,yrotModel);
        F = TriScatteredInterp(DT,sCalc);  sCalcGrid = F(X,Y);
        F = TriScatteredInterp(DT,bCalc);  bCalcGrid = F(X,Y);
        F = TriScatteredInterp(DT,uCalc);  uCalcGrid = F(X,Y);
        F = TriScatteredInterp(DT,vCalc);  vCalcGrid = F(X,Y);
        F = TriScatteredInterp(DT,wCalc);  wCalcGrid = F(X,Y);
        
        
        DT = DelaunayTri(xrotMeas,yrotMeas);
        F = TriScatteredInterp(DT,sData);  sDataGrid = F(X,Y);
        F = TriScatteredInterp(DT,bData);  bDataGrid = F(X,Y);
        F = TriScatteredInterp(DT,uData);  uDataGrid = F(X,Y);
        F = TriScatteredInterp(DT,vData);  vDataGrid = F(X,Y);
        F = TriScatteredInterp(DT,wData);  wDataGrid = F(X,Y);
        
        % replace NaNs with zeroth order values
        uDataGrid(isnan(uDataGrid))=U;
        vDataGrid(isnan(vDataGrid))=0;
        wDataGrid(isnan(wDataGrid))=0;
        sDataGrid(isnan(sDataGrid))=sMeasMean;
        bDataGrid(isnan(bDataGrid))=bMeasMean;
        
        uCalcGrid(isnan(uCalcGrid))=U;
        vCalcGrid(isnan(vCalcGrid))=0;
        wCalcGrid(isnan(wCalcGrid))=0;
        sCalcGrid(isnan(sCalcGrid))=sMeasMean;
        bCalcGrid(isnan(bCalcGrid))=bMeasMean;
        
        %
        %  figure ;surface(X,Y,sCalcGrid) ; title('sCalcGrid') ; xlabel('x') ; ylabel('y');
        %  figure ;surface(X,Y,bCalcGrid) ; title('bCalcGrid') ; xlabel('x') ; ylabel('y');
        %  figure ;surface(X,Y,uCalcGrid) ; title('uCalcGrid') ; xlabel('x') ; ylabel('y');
        %  figure ;surface(X,Y,vCalcGrid) ; title('vCalcGrid') ; xlabel('x') ; ylabel('y');
        %  figure ;surface(X,Y,wCalcGrid) ; title('wCalcGrid') ; xlabel('x') ; ylabel('y');
        %
        %  figure ;surface(X,Y,sDataGrid) ; title('sDataGrid') ; xlabel('x') ; ylabel('y');
        %  figure ;surface(X,Y,bDataGrid) ; title('bDataGrid') ; xlabel('x') ; ylabel('y');
        %  figure ;surface(X,Y,uDataGrid) ; title('uDataGrid') ; xlabel('x') ; ylabel('y');
        %  figure ;surface(X,Y,vDataGrid) ; title('vDataGrid') ; xlabel('x') ; ylabel('y');
        %  figure ;surface(X,Y,wDataGrid) ; title('wDataGrid') ; xlabel('x') ; ylabel('y');
        
        % I want to invert for differences:
        %
        du=uDataGrid-uCalcGrid; dv=vDataGrid-vCalcGrid; dw=wDataGrid-wCalcGrid; ds=sDataGrid-sCalcGrid;
        
        
        
        
        
        
        %
        %     estimate x with y given and y = K x
        %
        %     m measurements of each component, ie m=3 nx ny
        %     n elements of each component of the system state, ie n=2 nx ny
        %   [ds]=[Tsb  Tsc ]  [db]
        %   [du] [Tub  Tuc ] [dc]
        %   [dv] [Tvb  Tvc ]
        %   [dw] [Twb  Twc ]
        %
        %    Solution (n-form) ; x=xa+(K' iSe K + iSa )\ K' iSe (y-Kxa)
        %
        %  using Ce =[Cse 0  0  0 ]    diagonal matrices, ie Cse=seps^2 ones(
        %            [ 0 Cue 0  0 ]
        %            [ 0  0 Cve 0 ]
        %            [ 0  0  0 Cwe]
        %
        %   I need iSa in Fourier space:
        %   I use a pragmatic approach where iSa=L where L is the negative 2D discreate Laplacian
        %   this is using Sobolev H^1 penalty functional
        %
        %  iSa=[gb L    0]  where gb and gc are weighting parameters
        %      [0    gc L]
        %
        %
        % in the case of only one element in the system state, for exampel x=b we find
        % The 2D discreate neg. Laplacian : L b where L is (nx ny) x (nx ny) and b=vec(B) where B is the nx ny matrix
        % of bed elevation values
        % Using periodic boundary conditions L is bccb and is compleatly defined through its first column
        % and L b=ifft2(fft2(R) .* fft2(X)) where R is the nx x ny matrix obtained by the operaton
        % array(L.,1), ie be R=reshape(L(:,1),nx,ny); and X=array(x) , ie X=reshape(x,nx,ny) where x is the nx ny vector
        % of bed elevation values
        
        %[ testing
        if GenForwardData==1
            clear all
            H=1000;
            xmin=-100*H ; xmax=100*H; ymin=-100*H ; ymax=100*H;
            nx=64 ; ny=64;
            time=NaN;
            x=linspace(xmin,xmax,nx); y=linspace(ymin,ymax,ny);  dx=x(2)-x(1);dy=y(2)-y(1);
            
            %% get synthetic data
            [ds,du,dv,dw,db,dc,ds0] = SynthData2D(x,y,time);
            
            [X,Y]=ndgrid(x,y);
            rho=917; g=9.81/1000; alpha0=0.01;
            m=1 ; n=1 ; SlipRatio=1000; ud=1 ; taub=rho*g*H*sin(alpha0); AGlen=ud/(2*taub^n*H/(n+1));
            ub=SlipRatio*ud;
            C0=ub/taub^m;  eta=1/(2*AGlen);
            
            b=-H+db; c=C0*(1+dc);
            
            % add noise
            seps=0.0001*H; ueps=0.0001*mean(du(:)); veps=ueps; weps=ueps;
            seps=1e-8*H; ueps=1e-8*mean(du(:)); veps=ueps; weps=ueps;
            gs=1/seps^2; gu=1/ueps^2; gv=1/veps^2; gw=1/weps^2 ;
            ds=ds+seps*rand(nx,ny);
            du=du+ueps*rand(nx,ny);
            dv=dv+veps*rand(nx,ny);
            dw=dw+weps*rand(nx,ny);
        end
        
%         
%         figure ;contourf(X,Y,du-mean(du(:))) ; title('du') ; xlabel('x') ; ylabel('y'); colorbar
%         figure ;contourf(X,Y,dv) ; title('dv') ; xlabel('x') ; ylabel('y'); colorbar
%         figure ;contourf(X,Y,dw) ; title('dw') ; xlabel('x') ; ylabel('y'); colorbar
%         figure ;contourf(X,Y,ds) ; title('ds') ; xlabel('x') ; ylabel('y'); colorbar
        %figure ; quiver(X,Y,du-mean(du(:)),dv) ; title(' du dv')
        
        
        
        
        [Tsb,Tsc,Tss0,Tub,Tuc,Tus0,Tvb,Tvc,Tvs0,Twb,Twc,Tws0] = ForwardFunctions2D(dx,dy,nx,ny,time,alpha0,H,eta,C0,rho,g,m);
        
        %% Generate 2D Fourier representer of the 2D negative Laplacian with periodic boundary conditions
        lc=[4 -1 zeros(1,nx-3) -1 -1 zeros(1,nx-1) zeros(1,nx*(ny-3)) -1 zeros(1,nx-1)] ; %/(dx*dx+dy*dy)/2;
        R=ifft2(reshape(lc,nx,ny));
        
        
        
        %% Solves for b and c using s and u (vectorized)
        fftu=ifft2(du) ; ffts=ifft2(ds); fftv=ifft2(dv);  fftw=ifft2(dw);
        
        %  A11=gs.*conj(Tsb).*Tsb+gu.*conj(Tub).*Tub+gv.*conj(Tvb).*Tvb+gb.*R;
        %  A21=gs.*conj(Tsc).*Tsb+gu.*conj(Tuc).*Tub+gv.*conj(Tvc).*Tvb;
        %  A12=gs.*conj(Tsb).*Tsc+gu.*conj(Tub).*Tuc+gv.*conj(Tvb).*Tvc;
        %  A22=gs.*conj(Tsc).*Tsc+gu.*conj(Tuc).*Tuc+gv.*conj(Tvc).*Tvc+gc.*R;
        %  b1=gs.*conj(Tsb).*ffts+gu.*conj(Tub).*fftu+gv.*conj(Tvb).*fftv;
        %  b2=gs.*conj(Tsc).*ffts+gu.*conj(Tuc).*fftu+gv.*conj(Tvc).*fftv;
        
        
        A11=gs.*conj(Tsb).*Tsb+gu.*conj(Tub).*Tub+gv.*conj(Tvb).*Tvb+gw.*conj(Twb).*Twb+gb.*R;
        A21=gs.*conj(Tsc).*Tsb+gu.*conj(Tuc).*Tub+gv.*conj(Tvc).*Tvb+gw.*conj(Twc).*Twb;
        A12=gs.*conj(Tsb).*Tsc+gu.*conj(Tub).*Tuc+gv.*conj(Tvb).*Tvc+gw.*conj(Twb).*Twc;
        A22=gs.*conj(Tsc).*Tsc+gu.*conj(Tuc).*Tuc+gv.*conj(Tvc).*Tvc+gw.*conj(Twc).*Twc+gc.*R;
        b1=gs.*conj(Tsb).*ffts+gu.*conj(Tub).*fftu+gv.*conj(Tvb).*fftv+gw.*conj(Twb).*fftw;
        b2=gs.*conj(Tsc).*ffts+gu.*conj(Tuc).*fftu+gv.*conj(Tvc).*fftv+gw.*conj(Twc).*fftw;
        
        detA=A11.*A22-A21.*A12;
        invA11=A22./detA;
        invA12=-A12./detA;
        invA21=-A21./detA;
        invA22=A11./detA;
        db_est=invA11.*b1+invA12.*b2;
        dc_est=invA21.*b1+invA22.*b2;
        db_est=real(fft2(db_est)) ; dc_est=real(fft2(dc_est));
        db_est=db_est-mean(db_est(:)) ; dc_est=dc_est-mean(dc_est(:));
        
        
        %b_est=bMeasMean+db_est ;
        %c_est=C0*(1+dc_est);
        
        %
        %  figure ; contourf(X,Y,real(db_est)) ; title(' db est from s, u, v  ') ; colorbar
        %  figure ; contourf(X,Y,real(dc_est)) ; title(' dc est from s, u, v  ') ; colorbar
        %
        if GenForwardData==1
            subplot(2,2,1) ; contourf(X,Y,real(db)) ; title(' db  ') ; colorbar
            subplot(2,2,2) ;contourf(X,Y,real(dc)) ; title(' dc   ') ; colorbar
            subplot(2,2,3) ; contourf(X,Y,real(db_est-db)) ; title(' b_est-b ') ; colorbar
            subplot(2,2,4) ;contourf(X,Y,real(dc_est-dc)) ; title(' c_est-c  ') ; colorbar
        end
        
        % interpolate estimated values onto the numerical grid
        
        Xrot=X*cos(theta)-Y*sin(theta); Yrot=X*sin(theta)+Y*cos(theta);
        DT = DelaunayTri(Xrot(:),Yrot(:));
        F = TriScatteredInterp(DT,db_est(:));  db = F(xModel,yModel);
        F = TriScatteredInterp(DT,dc_est(:));  dc = F(xModel,yModel);
        
        TRI= delaunay(xModel,yModel);
%         figure ; trisurf(TRI,xModel,yModel,db) ;  title('db')
%         figure ; trisurf(TRI,xModel,yModel,dc) ;  title('dc')
%         
        
        % only update in the centre
        xc=0 ; yc=0 ; Radius=25000 ;
         distance=sqrt((xModel-xc).^2+(yModel-yc).^2);
        ind=find(distance < Radius);
        
        bModel(ind)=bModel(ind)+db(ind);
        C(ind)=C(ind).*(1+dc(ind));
        figure ; trisurf(TRI,xModel,yModel,bModel-bModelStart) ;  title('bModel-bModelStart')
        figure ; trisurf(TRI,xModel,yModel,C-CModelStart) ;  title('C-CModelStart')
        
    end
    
     save Deblurr2D time ...
        sModel uModel vModel bModel BModel ...
        sMeas uMeas vMeas wMeas bMeas BMeas xMeas yMeas ...
        Experiment ...
        coordinates connectivity Nnodes Nele nip nod etaInt gfint AGlen C ...
        Luv Luvrhs lambdau n m alpha rho rhow g Itime bModelStart CModelStart
    
    
    
end


