
function Imisfit=ModelDataMisfit(sModel,bModel,hModel,uModel,vModel,B,M,sMeas,uMeas,vMeas,wMeasInt,bMeas,BMeas,coordinates,connectivity,nip,Itype)
    % misfit calculated as integral (not vectorized)
    
    % for plotting misfit I need BModel, must add at some sate
    
    Nnodes=length(coordinates); %[Nele,nod]=size(connectivity);
    PlotMisfit=0; 
    
    if Itype==1  % discrete
        
        Imisfit=sum((uModel-uMeas).^2)+sum((vModel-vMeas).^2);
        
    elseif Itype==10
        
        %Imisfit=sum((uModel-uMeas).^2)+sum((vModel-vMeas).^2)+sum((wModelInt-wMeasInt).^2);
        k=2*Nnodes;
        B=[sparse(1:k,1:k,1,k,k) ; B];
        Bu=B*[uModel;vModel] ; d=[uMeas;vMeas;wMeasInt] ;
        Imisfit=(Bu-d)'*M*(Bu-d);
        
        
    else
        
        Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
        ndim=2; dof=2; neq=dof*Nnodes;
        
        [points,weights]=sample('triangle',nip,ndim);
        funInt=cell(nip); derInt=cell(nip);
        for Iint=1:nip
            funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
            derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
        end
        
        Imisfit=0;
        
        for Iele=1:Nele
            con=connectivity(Iele,:);  % nodes of element
            coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
            
            uModel_l=uModel(connectivity(Iele,:)); vModel_l=vModel(connectivity(Iele,:)) ;
            uMeas_l=uMeas(connectivity(Iele,:)); vMeas_l=vMeas(connectivity(Iele,:)) ;
            
            % g_l=connectivity(Iele,:);
            
            
            for Iint=1:nip                           % loop over integration points
                
                fun=funInt{Iint} ;
                der=derInt{Iint};
                J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
                detJ=det(J);  % det(dof x dof) matrix
                detJw=detJ*weights(Iint);
                
                uModelI=uModel_l'*fun ; vModelI=vModel_l'*fun ;
                uMeasI=uMeas_l'*fun ; vMeasI=vMeas_l'*fun ;
                F=(uModelI-uMeasI)^2+(vModelI-vMeasI)^2;
                Imisfit=Imisfit+F*detJw ;
                
            end
        end
        
    end
    
    Imisfit=real(Imisfit);
    
    if PlotMisfit==1
        x=coordinates(:,1); y=coordinates(:,2); DTxy = DelaunayTri(x,y); TRIxy=DTxy.Triangulation;
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

