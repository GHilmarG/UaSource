function [db,dc] = Adjoint2D(NaN,BoundaryNodes,...
        sModel,uModel,vModel,hModel,BModel,b_prior,...
        sMeas,uMeas,vMeas,wMeasInt,bMeas,BMeas,xMeas,yMeas,...
        Experiment,...
        coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,AGlen,C,C_prior,...
        Luv,Luvrhs,lambdau,n,m,alpha,rho,rhow,g,Itime);
    
    global solutionphase
    
   % does not converge to correct solutision
   % check what the residuals are doing
   % how does including the  vertical velocity in the misfit affect the solution
    
    restart=0; icountmax=100; nsteepestdescent=25; 
    PlotMisfit=1;
    conjgrad=1;      % 1 if conjugated gradients are used 
    TestDirection=0; % tests if the cost function really does decrease in the direction of neg. gradient by calculating a few values 
    
    disp(' Adjoint2D')

    x=coordinates(:,1); y=coordinates(:,2); DTxy = DelaunayTri(x,y); TRIxy=DTxy.Triangulation;
    
    nInt=numel(etaInt); % number of integration points
    %figure(900) ; trisurf(TRIxy,x,y,hModel-(sMeas-bMeas)) ;  title(' hModel-(sMeas-bMeas)')
    %figure(901) ; trisurf(TRIxy,x,y,sModel-sMeas) ;  title(' sModel-sMeas')
     
     
    icount=1;
    % J  is the objective function,
    % J= Imisfit  + gamma Ireg
    
    iModelType=1;  % discrete norm
    iMisfitType=10; % discrete, 1 without and 10 with w as well
    
    
    beta= sqrt(C.^(-1/m).* (sqrt(uModel.*uModel+vModel.*vModel)).^(1/m-1));
    beta_prior= sqrt(C_prior.^(-1/m).* (sqrt(uModel.*uModel+vModel.*vModel)).^(1/m-1));
    beta_error=mean(beta);
    beta_error=1.0;
    
    [b_prior]=FunctionBedrock2d('forwardC',coordinates); % b_prior is const. sloping bed
    b_error=300;
    
    nM=2*Nnodes+nInt;
    
    uError=zeros(Nnodes,1)+1; vError=zeros(Nnodes,1)+1; wError=zeros(nInt,1)+1;
    M=sparse(1:nM,1:nM,[1./uError.^2;1./vError.^2;1./wError.^2],nM,nM);   % I=(B u-d)' M (B u -d ) 
        
    n1beta=1 ; n2beta=length(beta) ; n1b=n2beta+1 ; n2b=n1b+length(sModel)-1;
    
    if restart==1
        disp(' Adjoint Restart ')
        load AdjointRestart q Jvector AngleVector l ;
        Jlength=length(Jvector) ;  Jvector=[Jvector; zeros(icountmax,1)];
        Anglelength=length(AngleVector) ;  AngleVector=[AngleVector; zeros(icountmax,1)];
        beta=q(n1beta:n2beta);
        bModel=q(n1b:n2b);
        hModel=sModel-bModel;
    else
        Jlength=1; Jvector=zeros(icountmax,1); 
        Anglelength=1; AngleVector=zeros(icountmax,1); 
        l=1;
        bModel=sModel-hModel;
        q=[beta;bModel] ; % q ist the distributed parameter to be estimated
    end
    % iType=0;  % square of difference between beta and beta_prior
    
    disp(' Calculating cost function: includes a LS misfit term and  regularisation terms ')
    
   

    
    
    fJpar={beta_prior,beta_error,b_prior,b_error,M,iModelType,iMisfitType, ...
        sModel,uModel,vModel,BModel,sMeas,uMeas,vMeas,wMeasInt,bMeas,BMeas,...
        coordinates,connectivity,nip,etaInt,gfint,AGlen,Luv,Luvrhs,n,m,alpha,rho,rhow,g,icount,n1beta,n2beta,n1b,n2b};
    
     fdJdqPar={beta_prior,beta_error,b_prior,b_error,M,iModelType,iMisfitType,...
        BoundaryNodes,uModel,vModel,sModel,BModel,...
        uMeas,vMeas,wMeasInt,coordinates,connectivity,nip,...
        etaInt,gfint,AGlen,Luv,Luvrhs,n,m,alpha,rho,rhow,g,icount,n1beta,n2beta,n1b,n2b};
    
    disp('Initial value of the cost function J ')
    J=fJ(q,fJpar{:});
    
    if restart~=1 ;  Jvector(1)=J; end
    
    disp('Initial search direction: dJ/dq')
     
    dJdq=fdJdq(q,fdJdqPar{:});


  
   
    direction=-dJdq;
    icount=1; ftol=1e-5;
    

    
    
    if TestDirection==1
        % check if direction is sensible
        
        figure(201) ; trisurf(TRIxy,x,y,hModel) ;  title(' hModel')
        figure(202) ; trisurf(TRIxy,x,y,bModel) ;  title(' bModel')
        figure(200) ; trisurf(TRIxy,x,y,b_prior) ;  title(' b prior')
        figure(203) ; trisurf(TRIxy,x,y,bModel-b_prior) ;  title('bModel- b prior')
        dl=0.5; NN=20;
        for I=1:NN
            l=(I-1-NN/4)*dl;
            J=f1dim(l,@(dummy) fJ(dummy,fJpar{:}), q,direction);
            disp([' J : ',num2str(J),', l : ',num2str(l)])
            Jvec(I)=J ; lvec(I)=l;
        end
        
        figure(999) ; plot(lvec,Jvec,'-o');

    end


   
    
    
    while icount<= icountmax
        
        ax=0 ; bx=l/10;   Jlast=J;
        
     
        
        fJpar={beta_prior,beta_error,b_prior,b_error,M,iModelType,iMisfitType, ...
            sModel,uModel,vModel,BModel,sMeas,uMeas,vMeas,wMeasInt,bMeas,BMeas,...
            coordinates,connectivity,nip,etaInt,gfint,AGlen,Luv,Luvrhs,n,m,alpha,rho,rhow,g,icount,n1beta,n2beta,n1b,n2b};
        
        fdJdqPar={beta_prior,beta_error,b_prior,b_error,M,iModelType,iMisfitType,...
            BoundaryNodes,uModel,vModel,sModel,BModel,...
            uMeas,vMeas,wMeasInt,coordinates,connectivity,nip,...
            etaInt,gfint,AGlen,Luv,Luvrhs,n,m,alpha,rho,rhow,g,icount,n1beta,n2beta,n1b,n2b};
        
                
        % function  [ax,bx,cx,fa,fb,fc]=mnbrak(func,ax,bx)
        %function [ result ] = f1dim(l,func,p,n)
        
        [ax,bx,cx,fa,fb,fc]=mnbrak(@(l) f1dim(l,@(dummy) fJ(dummy,fJpar{:}), q,direction), ax,bx);

        disp([' ax = ',num2str(ax),', bx = ',num2str(bx),', cx = ',num2str(cx)])
        disp([' fa = ',num2str(fa),', fb = ',num2str(fb),', fc = ',num2str(fc)])
        hold off ;figure(10) ; plot([ax bx cx ],[fa fb fc]) ; title(' mnbrak ') ; xlabel(' x '); ylabel('func') ; hold on
        
        
        % so the best value after mnbrack is q=q+bx
%         qTest=q+bx*direction;
%         betaTest=qTest(n1beta:n2beta);
%         bModelTest=qTest(n1b:n2b);
%         figure(201) ; trisurf(TRIxy,x,y,betaTest) ;  title(' betaTest')
%         figure(202) ; trisurf(TRIxy,x,y,bModelTest-bModel) ;  title(' bModelTest-bModel')
%         figure(203) ; trisurf(TRIxy,x,y,bModel) ;  title(' bModel')
        
     
        
        %function [funmin,xmin]=dbrentghg(func,dfunc,ax,bx,cx,fa,fb,fc,tol,funlast)
        %function [ result ] = f1dim(l,func,p,n)
        disp(' dbrentghg ' )
        [J,l] = dbrentghg(@(dummy1) f1dim(dummy1,@(dummy2) fJ(dummy2,fJpar{:}),q,direction),...
            @(dummy3) df1dim(dummy3,@(dummy4) fdJdq(dummy4,fdJdqPar{:}),...
            q,direction),...
            ax,bx,cx,fa,fb,fc,1e-3,Jlast);
        
        
        
        Jvector(Jlength+icount)=J;
        disp([' J before = ',num2str(Jlast),' after = ',num2str(J),', before-after = ',num2str(Jlast-J),' Jafter/Jbefore : ',num2str(J/Jlast)])
        
        hold on ; figure(10) ; plot(l,J,'+') ; hold off
        
        q=q+bx*direction;
        
        bModelLast=bModel; betaLast=beta;
        beta=q(n1beta:n2beta);  bModel=q(n1b:n2b);
      
        figure(210) ; trisurf(TRIxy,x,y,beta) ;  title(' beta')
        figure(220) ; trisurf(TRIxy,x,y,bModel) ;  title(' bModel')
        figure(2000) ; trisurf(TRIxy,x,y,bModel-b_prior) ;  title(' bModel-b prior')
        figure(230) ; trisurf(TRIxy,x,y,bModel-bModelLast) ;  title(' bModel-bModelLast')
        
        
        figure(3000) ; trisurf(TRIxy,x,y,beta-beta_prior) ;  title(' beta-beta prior')
        figure(240) ; trisurf(TRIxy,x,y,beta-betaLast) ;  title(' beta-betaLast')
        
        
%         if 2*abs(J-Jlast) < ftol*(abs(J)+abs(Jlast))+100*eps ;
%             disp([' non-linear cg converged after iteration # ', num2str(icount)]) ;
%             break ;
%         end
%         
        
        
        dJdqlast=dJdq;
        dJdq=fdJdq(q,fdJdqPar{:});
         
                
        
        % new conjucated gradient direction
        if conjgrad==1 && (icount > 1 || restart==1 )
            
            teta=dJdq'*dJdq/(dJdqlast'*dJdqlast);   % Fletcher and Reeves
            
            
            %teta=(dJdq-dJdqlast)'*dJdq/(dJdqlast'*dJdqlast);  % Polak and Ribiere
            
             if mod(icount,nsteepestdescent)==0 ;
                teta=0 ;% apparantly it is good to operate cg in cycles (or so I'm told)
                disp(' resetting cj ')
             end

            
             direction=-dJdq + teta * direction ; % new search direction
             
             % angle between steepest decent and cg
             angle=180*acos(-dJdq'*direction/norm(dJdq)/norm(direction))/pi;
             
             AngleVector(Anglelength+icount)=angle;
             
             disp([' teta : ', num2str(teta),' angle ' ,num2str(angle)])
             if teta < 0 ; disp(' note teta < 0 !!!!') ; end
             
             
            
        else  % do the first step using steepest decent
            direction=-dJdq;
        end
        icount=icount+1;
        
        save AdjointRestart q Jvector AngleVector l ;
        
    end
    
    
    
    if PlotMisfit==1
        
        
        C= ((sqrt(uModel.*uModel+vModel.*vModel)).^(1-m)).*beta.^(-2) ;
        
        % disp(' Step#1: Solve the forward model ')
        
        
        
        solutionphase=1 ; % solving state equation
        [uModel,vModel,lambdau,kv]=SSTREAM2d(sModel,hModel,uModel,vModel,coordinates,connectivity,nip,...
            etaInt,gfint,AGlen,C,Luv,Luvrhs,lambdau,n,m,alpha,rho,rhow,g,icount);
        
        
        
        figure(701) ; trisurf(TRIxy,x,y,uModel) ;  title(' uModel')
        figure(702) ; trisurf(TRIxy,x,y,uMeas) ;  title(' uMeas')
        figure(703) ; trisurf(TRIxy,x,y,uModel-uMeas) ;  title(' uModel-uMeas')
        
        figure(711) ; trisurf(TRIxy,x,y,vModel) ;  title(' vModel')
        figure(712) ; trisurf(TRIxy,x,y,vMeas) ;  title(' vMeas')
        figure(713) ; trisurf(TRIxy,x,y,vModel-vMeas) ;  title(' vModel-vMeas')
     
        figure(721) ; trisurf(TRIxy,x,y,sModel) ;  title(' sModel')
        figure(722) ; trisurf(TRIxy,x,y,sMeas) ;  title(' sMeas')
        figure(723) ; trisurf(TRIxy,x,y,sModel-sMeas) ;  title(' sModel-sMeas')
        
        figure(731) ; trisurf(TRIxy,x,y,bModel) ;  title(' bModel')
        figure(732) ; trisurf(TRIxy,x,y,bMeas) ;  title(' bMeas')
        figure(733) ; trisurf(TRIxy,x,y,bModel-bMeas) ;  title(' bModel-bMeas')
        
        hMeas=sMeas-bMeas;
        figure(741) ; trisurf(TRIxy,x,y,hModel) ;  title(' hModel')
        figure(742) ; trisurf(TRIxy,x,y,hMeas) ;  title(' hMeas')
        figure(743) ; trisurf(TRIxy,x,y,hModel-hMeas) ;  title(' hModel-hMeas')
        
        
        
        [wModelInt,B,xint,yint] = VertVelMatrixVector(hModel,bModel,uModel,vModel,coordinates,connectivity,nip);
         
        DTint = DelaunayTri(xint,yint); TRIint=DTint.Triangulation;
        figure(751) ; trisurf(TRIint,xint,yint,wModelInt) ;  title(' wModelInt')
        figure(752) ; trisurf(TRIint,xint,yint,wMeasInt) ;  title(' wMeasInt')
        figure(753) ; trisurf(TRIint,xint,yint,wModelInt-wMeasInt) ;  title(' wModelInt-wMeasInt')
        
        
    end


    
    
    
    Cretrieved= ((sqrt(uModel.*uModel+vModel.*vModel)).^(1-m)).*beta.^(-2) ;
    [Ctrue]=FunctionSlipperiness2d('forwardb',coordinates);
    figure(103) ; trisurf(TRIxy,x,y,Cretrieved) ;  title(' C retrieved ')
    figure(104) ; trisurf(TRIxy,x,y,Ctrue) ;  title(' C true ')
    figure(105) ; trisurf(TRIxy,x,y,Ctrue-Cretrieved) ;  title(' C true - C retrieved ')
    
    figure(50) ; semilogy(Jvector,'-o')
    figure(51) ; plot(Jvector,'-o')
    figure(52) ; plot(AngleVector,'-o') ; title('angle between steepest decent and gc ')
    
    
    
  
    
    
    error('safsad')
end

