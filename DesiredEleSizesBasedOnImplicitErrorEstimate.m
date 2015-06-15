
function  [coordinates,connectivity]=DesiredEleSizesBasedOnImplicitErrorEstimate(Experiment,MeshBoundaryCoordinates,Boundary0,...
        s0,b0,S0,B0,h0,u0,v0,coordinates0,connectivity0,nip,AGlen0,C0,Luv0,Luvrhs0,lambdauv0,n,m,alpha,rho,rhow,g,Itime,CtrlVar)
    
    
    infovector=zeros(20,5)+NaN; isScaledDown=0;
    
    fprintf(CtrlVar.fidlog,' DesiredEleSizesBasedOnImplicitErrorEstimate \n');
    
    MeshSizeMax=CtrlVar.MeshSizeMax ; MeshSizeMin=CtrlVar.MeshSizeMin ;
    
    %% if u and v has not been calculated on orignal grid, do so now
    
    [u0,v0,lambda0]= SSTREAM2dNR(s0,S0,B0,h0,u0,v0,coordinates0,connectivity0,Boundary0,nip,AGlen0,C0,Luv0,Luvrhs0,lambdauv0,n,m,alpha,rho,rhow,g,Itime,CtrlVar);
        
    
    x0=coordinates0(:,1) ; y0=coordinates0(:,2);
    [Nele0,nod0]=size(connectivity0);  Nnodes0=max(connectivity0(:));
    fprintf(CtrlVar.fidlog,' Meshrefinement on input \t # Elements: %-i \t \t # Nodes=%-i \n',Nele0,Nnodes0);
    % too do: after a certain number of local meshrefinements, estimate element sizes, smooth it and do a global mesh refinement
    
    NEleMax=CtrlVar.MaxNumberOfElements;
    
    iRefinement=0; JobVar.iMeshRefinement=0;
    while iRefinement < CtrlVar.AdaptMeshIterations && isScaledDown==0
        
        iRefinement=iRefinement+1;
        fprintf(CtrlVar.fidlog,'\n \n');
        size(coordinates0)
        [EleError0,NormalizedEleError0,rmsEleError0]=fe2dErrorEstimate(iRefinement,Experiment,coordinates0,connectivity0,...
            Boundary0,MeshBoundaryCoordinates,u0,v0,h0,s0,b0,AGlen0,C0,n,m,rho,rhow,g,Itime,nip,DTxy0,TRIxy0,...
            CtrlVar);  % broken due to external changes, presumably easy to fix
        
        % Needed here: a decision on whether to do element refinement/remeshing at all
        
        
        if CtrlVar.MeshRefinement==1
            
            NormalizedErrorCrit=0.01;
            EleRefine0=NormalizedEleError0> NormalizedErrorCrit;
            %EleRefine0=sum(RefineNodes0(connectivity0),2) > 1 ; % if more than one node fullfill criteria, subdivide element
            NEleRefine0=sum(EleRefine0);
            
            eFrac=0.1;  % max fraction of elements refined
            
            %if NEleRefine0/Nele0>eFracMax || NEleRefine0/Nele0<eFrac
            %   fprintf(CtrlVar.fidlog,' Refinement criteria changed so that no more than %-f percent of elements subdivided \n',100*eFrac)
            NormalizedError0Sorted=sort(NormalizedEleError0) ;
            NormalizedErrorCrit=NormalizedError0Sorted(floor(numel(NormalizedError0Sorted)*(1-eFrac))); % about 10% of elements refined at each step
            EleRefine0=NormalizedEleError0> NormalizedErrorCrit;
            NEleRefine0=sum(double(EleRefine0));
            %end
            
            fprintf(CtrlVar.fidlog,' Number of elements to be subdivided : %-i or %-f percent of total number of elements \n',NEleRefine0,100*NEleRefine0/Nele0);
            
            
            if sum(EleRefine0)==0
                fprintf(CtrlVar.fidlog,' No elements to subdivide! \n');
                break
            end
            
            
            [coordinates0,connectivity0] = refine(coordinates0,connectivity0,EleRefine0);
            if  mod(JobVar.iMeshRefinement,1)==0
                fprintf(CtrlVar.fidlog,' smoothing mesh \n ');
                MaxSmoothIterations=1;
                [coordinates0,connectivity0] = smoothmesh(coordinates0,connectivity0,MaxSmoothIterations);
            end
            JobVar.iMeshRefinement=JobVar.iMeshRefinement+1;
            
            NnodesLast=Nnodes0 ; NeleLast=Nele0;
            Nnodes0=max(connectivity0(:)); [Nele0,nod0]=size(connectivity0);
            fprintf(CtrlVar.fidlog,' Local meshrefinement of original mesh %-2.2i \t # Elements: %-i \t \t # Nodes=%-i \t from # Ele: %-i \t Nodes: %-i\n',iRefinement,Nele0,Nnodes0,NeleLast,NnodesLast);
            
            % It is not good to do repeated mesh refinements
            
            if Nele0>NEleMax
                fprintf(CtrlVar.fidlog,' local mesh refinement has resulted Nele=%-i > NEleMax=%-i and is stopped \n',Nele0,NEleMax);
                iRefinement = CtrlVar.AdaptMeshIterations ;
            end
            %
            %         if mod(JobVar.iMeshRefinement,10)==0
            %             fprintf(CtrlVar.fidlog,' global meshing \n ')
            %             [EleArea0,xEleCentre0,yEleCentre0]=TriAreaFE(coordinates0,connectivity0);
            %             CtrlVar.MeshSize=zeros(numel(xEleCentre0(1:1:end)),3);
            %
            %             CtrlVar.MeshSize(:,1)=xEleCentre0(1:1:end) ; CtrlVar.MeshSize(:,2)=yEleCentre0(1:1:end); CtrlVar.MeshSize(:,3)=sqrt(2)*sqrt(EleArea0(1:1:end));
            %
            %
            %             if Nele0>NEleMax
            %                 CtrlVar.MeshSize(:,3)=CtrlVar.MeshSize(:,3)*sqrt(Nele0/NEleMax);
            %                 fprintf(CtrlVar.fidlog,'Nele0=%-i > NEleMax=%-i , hence desired meshsize is scaled down \n',Nele0,NEleMax)
            %             end
            %
            %             %DTele0 = DelaunayTri(xEleCentre0,yEleCentre0);
            %             %trisurf(DTele0.Triangulation,xEleCentre0,yEleCentre0,EleArea0) ;  title(' desired ele sizes ')
            %             [coordinates0,connectivity0]=genmesh2d(Experiment,MeshBoundaryCoordinates,CtrlVar);
            %             Nnodes0=max(connectivity0(:));
            %             [Nele0,~]=size(connectivity0);
            %
            %
            %             fprintf(CtrlVar.fidlog,' global remeshing %-2.2i \t # Elements: %-i \t \t # Nodes=%-i \n',iRefinement,Nele0,Nnodes0)
            %
            %             JobVar.iMeshRefinement=0;
            %         end
            
        else  % global remeshing
            fprintf(CtrlVar.fidlog,' global meshing \n ');
            
            [EleArea0,xEleCentre0,yEleCentre0]=TriAreaFE(coordinates0,connectivity0);
            CtrlVar.MeshSize=zeros(numel(xEleCentre0(1:1:end)),3);
            
            % err0= c d0^k
            % err1= c d1^k
            % -> err1=err0 (d1/d0)^k
            % -> d1= (err1/err0)^(1/k)   h0
            % mesh is equalibrated if err1 is same for all elements
            % hence d1= K (err1/err0)^(1/k) d0, where err1 is the desired elementwise error
            % and err0 is the estimated error for each element (which will vary from element to element as long as the
            % mesh is not equalibrated)
            %
            % err1 can be set if different ways for example: err1=gamma mean(err0)
            % where gamma is the desired overall reduction in error at each step
            % repeated until mean(err0) less than a given error, or elements too large, etc
            
            
            switch nod0
                case 3
                    k=1;
                case 6
                    k=2;
                case 10
                    k=3;
            end
            
            d0=sqrt(2*EleArea0);
            gamma=1; % ratio between new mean of new and old element error
            % if set to one, then error is simply equlibriated over the new mesh
            fEle=2;  % don't allow elements to be reduced by more that this factor
            
            % infovector
            
            EleError1=gamma*mean(EleError0);
            d1=d0.*(EleError1./EleError0).^(1/k);
            d1=max(d1,d0/fEle);
            EleSize=d1;
            
            
            EleSize=max(MeshSizeMin,EleSize);  % make sure that elements are not smaller than MeshSizeMin
            EleSize=min(MeshSizeMax,EleSize);  % and not larger than MeshSizeMax
            
            CtrlVar.MeshSize(:,1)=xEleCentre0(1:1:end) ; CtrlVar.MeshSize(:,2)=yEleCentre0(1:1:end);
            CtrlVar.MeshSize(:,3)=EleSize;
            %DTele0 = DelaunayTri(xEleCentre0,yEleCentre0);
            %figure ; trisurf(DTele0.Triangulation,xEleCentre0,yEleCentre0,EleArea0) ;  title(' ele dimension 0 ')
            %figure ; trisurf(DTele0.Triangulation,xEleCentre0,yEleCentre0,EleSize) ;  title(' desired new ele dimentions ')
            
            
            
            [coordinates0,connectivity0]=genmesh2d(Experiment,MeshBoundaryCoordinates,CtrlVar);
            NnodesLast=Nnodes0 ; NeleLast=Nele0;
            Nnodes0=max(connectivity0(:));  [Nele0,~]=size(connectivity0);
            
            
            while Nele0>1.2*NEleMax
                CtrlVar.MeshSize(:,3)=CtrlVar.MeshSize(:,3)*sqrt(Nele0/NEleMax);
                fprintf(CtrlVar.fidlog,'Nele0=%-i > NEleMax=%-i , hence desired meshsize is scaled down \n',Nele0,NEleMax);
                [coordinates0,connectivity0]=genmesh2d(Experiment,MeshBoundaryCoordinates,CtrlVar);
                Nnodes0=max(connectivity0(:));  [Nele0,~]=size(connectivity0);
                fprintf(CtrlVar.fidlog,'new Nele0 after scale down is %-i  \n',Nele0);
                isScaledDown=1;
            end
            
            
            fprintf(CtrlVar.fidlog,' Global meshrefinement of original mesh %-2.2i \t # Elements: %-i \t \t # Nodes=%-i \t from # Ele: %-i \t Nodes: %-i\n',iRefinement,Nele0,Nnodes0,NeleLast,NnodesLast);
            JobVar.iMeshRefinement=0;
            
            infovector(iRefinement,1)=rmsEleError0 ;
            infovector(iRefinement,2)=mean(EleError0) ;
            infovector(iRefinement,3)=std(EleError0) ;
            infovector(iRefinement,4)=NeleLast ;
            infovector(iRefinement,5)=Nele0 ;
            
        end
        
%         [DTxy0,TRIxy0,u0,v0,h0,s0,b0,S0,B0,AGlen0,C0,rho,Boundary0,...
%             ufixednode0,ufixedvalue0,vfixednode0,vfixedvalue0,utiedA0,utiedB0,vtiedA0,vtiedB0,hfixednode0,hfixedvalue0,htiedA0,htiedB0,...
%             Luv0,Luvrhs0,lambdauv0,Lh0,Lhrhs0,lambdah0]=...
%             DiagnosticSolutionGrid1toGrid2(Experiment,coordinates0,connectivity0,...
%             Boundary0,MeshBoundaryCoordinates,u0,v0,h0,s0,b0,AGlen0,C0,n,m,rho,rhow,g,Itime,nip,DTxy0,...
%             CtrlVar);
        
        x0=coordinates0(:,1) ; y0=coordinates0(:,2);
        
    end
    
    
    % update all variables
    coordinates=coordinates0 ; connectivity=connectivity0;
    
%     u=u0 ; v=v0 ; h=h0 ; s=s0 ;b=b0 ; S=S0 ; B=B0 ;AGlen=AGlen0 ; C=C0 ; Boundary=Boundary0 ;
%    
%     
%     ufixednode=ufixednode0 ; ufixedvalue=ufixedvalue0 ; vfixednode=vfixednode0 ; vfixedvalue= vfixedvalue0 ; utiedA=utiedA0 ; utiedB=utiedB0 ; vtiedA=vtiedA0 ; vtiedB=vtiedB0;
%     hfixednode=hfixednode0 ; hfixedvalue=hfixedvalue0 ; htiedA=htiedA0 ; htiedB=htiedB0;
%     Luv=Luv0; Luvrhs=Luvrhs0 ; lambdauv=lambdauv0; Lh=Lh0 ; Lhrhs=Lhrhs0 ; lambdah=lambdah0 ;
    
    
    %save TestSave ufixedvalue vfixedvalue ufixednode vfixednode u v coordinates connectivity DTxy TRIxy
    
    %save DesiredEleSizesBasedOnImplicitErrorEstimateInfo infovector
    
end


