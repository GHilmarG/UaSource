function [R,K,T,F]=KRTFuvhGeneralTG3...
        (CtrlVar,MUA,u,v,h,S,B,u0,v0,h0,as0,ab0,as1,ab1,dudt,dvdt,dt,AGlen,n,C,m,alpha,rho,rhow,g)
        
    % uvhAssembly
    if nargin~=25
        error('wrong number of input arguments ')
    end
        
    
    %
    % K= [Kxu Kxv Kxh]
    %    [Kyu Kyv Kyh]
    %    [Khu Khv Khh]
   
    persistent nSave
    
    
    if nargout==1 
        Ronly=1; 
    else
        Ronly=0;
    end

    if any(isnan(u)) ;  fprintf(CtrlVar.fidlog,' NaN in u on input to KRTFuvhGeneralTG3 \n'); end
    if any(isnan(v)) ;  fprintf(CtrlVar.fidlog,' NaN in v on input to KRTFuvhGeneralTG3 \n'); end
    if any(isnan(h)) ;  fprintf(CtrlVar.fidlog,' NaN in h on input to KRTFuvhGeneralTG3 \n'); end
    if any(isnan(u0)) ;  fprintf(CtrlVar.fidlog,' NaN in u0 on input to KRTFuvhGeneralTG3 \n'); end
    if any(isnan(v0)) ;  fprintf(CtrlVar.fidlog,' NaN in v0 on input to KRTFuvhGeneralTG3 \n'); end
    if any(isnan(h0)) ;  fprintf(CtrlVar.fidlog,' NaN in h0 on input to KRTFuvhGeneralTG3 \n'); end
    
    temp=CtrlVar.ResetThicknessToMinThickness;
    if ~CtrlVar.ResetThicknessInNonLinLoop
        CtrlVar.ResetThicknessToMinThickness=0;
    end
    
    [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
    
    if CtrlVar.MassBalanceGeometryFeedback>=2
        GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
        rdamp=CtrlVar.MassBalanceGeometryFeedbackDamping;
        if rdamp~=0
            as1Old=as1 ; ab1Old=ab1;
        end
        
        switch CtrlVar.MassBalanceGeometryFeedback
            case 2
                [as1,ab1]=DefineMassBalance(CtrlVar.UserVar,CtrlVar,MUA,CtrlVar.time+dt,s,b,h,S,B,rho,rhow,GF);
                dadh=zeros(MUA.Nnodes,1);
            case 3
                [as1,ab1,dasdh,dabdh]=DefineMassBalance(CtrlVar.UserVar,CtrlVar,MUA,CtrlVar.time+dt,s,b,h,S,B,rho,rhow,GF);
                dadh=dasdh+dabdh;
        end
        
        
        if rdamp~=0
            % I don't account for a potential dependency of as and ab
            % on h in the Hessian, so may need to dampen these changes
            as1=(1-rdamp)*as1+rdamp*as1Old;
            ab1=(1-rdamp)*ab1+rdamp*ab1Old;
        end
    else
        dadh=zeros(MUA.Nnodes,1);
    end
    
    
    
    
    CtrlVar.ResetThicknessToMinThickness=temp;
    
    [etaInt,~,~,exx,eyy,exy,Eint]=calcStrainRatesEtaInt(CtrlVar,MUA,u,v,AGlen,n);
    
    
   
    
    ndim=2;  neq=3*MUA.Nnodes;
    neqx=MUA.Nnodes ;

    if CtrlVar.TG3
        if isempty(nSave)
            nSave=1; fprintf(CtrlVar.fidlog,' Fully implicit uvh assembly based on third-order Taylor Galerkin \n');
        end
    end
  
    
    
    hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
    unod=reshape(u(MUA.connectivity,1),MUA.Nele,MUA.nod);
    vnod=reshape(v(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    if ~CtrlVar.CisElementBased 
        C=reshape(C(MUA.connectivity,1),MUA.Nele,MUA.nod);  % This could be called Cnod
    end
    
    
    
    Snod=reshape(S(MUA.connectivity,1),MUA.Nele,MUA.nod);
    Bnod=reshape(B(MUA.connectivity,1),MUA.Nele,MUA.nod);
    h0nod=reshape(h0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    u0nod=reshape(u0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    v0nod=reshape(v0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    as0nod=reshape(as0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    ab0nod=reshape(ab0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    as1nod=reshape(as1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    ab1nod=reshape(ab1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    dadhnod=reshape(dadh(MUA.connectivity,1),MUA.Nele,MUA.nod);
    rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
    bnod=reshape(b(MUA.connectivity,1),MUA.Nele,MUA.nod);
    dudtnod=reshape(dudt(MUA.connectivity,1),MUA.Nele,MUA.nod);
    dvdtnod=reshape(dvdt(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    
    
    ca=cos(alpha); sa=sin(alpha);
    
    
   % [points,weights]=sample('triangle',nip,ndim);
    
    if ~Ronly
        Kxu=zeros(MUA.Nele,MUA.nod,MUA.nod); Kxv=zeros(MUA.Nele,MUA.nod,MUA.nod);  Kxh=zeros(MUA.Nele,MUA.nod,MUA.nod);
        Kyu=zeros(MUA.Nele,MUA.nod,MUA.nod); Kyv=zeros(MUA.Nele,MUA.nod,MUA.nod);  Kyh=zeros(MUA.Nele,MUA.nod,MUA.nod);
        Khu=zeros(MUA.Nele,MUA.nod,MUA.nod); Khv=zeros(MUA.Nele,MUA.nod,MUA.nod);  Khh=zeros(MUA.Nele,MUA.nod,MUA.nod);
    else
        Kxu=[]; Kxv=[];  Kxh=[];
        Kyu=[]; Kyv=[];  Kyh=[];
        Khu=[]; Khv=[];  Khh=[];
    end
    
    Tx=zeros(MUA.Nele,MUA.nod);  Ty=zeros(MUA.Nele,MUA.nod); Fx=zeros(MUA.Nele,MUA.nod);  Fy=zeros(MUA.Nele,MUA.nod); Th=zeros(MUA.Nele,MUA.nod);  Fh=zeros(MUA.Nele,MUA.nod);
    
    Tx0=Tx ;Fx0=Fx; Ty0=Ty ; Fy0=Fy; Th0=Th ;Fh0=Fh;
    Kxu0=Kxu ; Kxv0=Kxv ; Kyu0=Kyu ; Kyv0=Kyv ; Kxh0=Kxh ; Kyh0=Kyh ; Khu0=Khu ; Khv0=Khv ; Khh0=Khh; 
    
   
    
    
    
    parfor Iint=1:MUA.nip
        
        [Tx1,Fx1,Ty1,Fy1,Th1,Fh1,Kxu1,Kxv1,Kyu1,Kyv1,Kxh1,Kyh1,Khu1,Khv1,Khh1]=...
            uvhAssembleIntPointImplicitTG3(Iint,ndim,MUA,...
            bnod,hnod,unod,vnod,C,h0nod,u0nod,v0nod,as0nod,ab0nod,as1nod,ab1nod,dudtnod,dvdtnod,dadhnod,Bnod,Snod,rhonod,...
            CtrlVar,rhow,g,m,etaInt,exx,eyy,exy,Eint,Ronly,ca,sa,dt,...
            Tx0,Fx0,Ty0,Fy0,Th0,Fh0,Kxu0,Kxv0,Kyu0,Kyv0,Kxh0,Kyh0,Khu0,Khv0,Khh0);
                 
        
        
        Tx=Tx+Tx1;  Fx=Fx+Fx1; 
        Ty=Ty+Ty1;  Fy=Fy+Fy1;
        Th=Th+Th1;  Fh=Fh+Fh1;
        
        Kxu=Kxu+Kxu1;        Kxv=Kxv+Kxv1;        
        Kyu=Kyu+Kyu1;        Kyv=Kyv+Kyv1;        
        Kxh=Kxh+Kxh1;        Kyh=Kyh+Kyh1;
        Khu=Khu+Khu1;        Khv=Khv+Khv1;        Khh=Khh+Khh1;
        
    end
    %%
    
%     if CtrlVar.InfoLevelCPU ;  
%         UaInfo.CPUuvhAssembly=UaInfo.CPUuvhAssembly+toc(tAssembly) ;  
%         UaInfo.CPUuvhAssemblyCounter=UaInfo.CPUuvhAssemblyCounter+1 ; 
%     end
%     
    
    %% assemble right-hand side
    
    T=sparseUA(neq,1); F=sparseUA(neq,1);
    
    for Inod=1:MUA.nod
        
        
        T=T+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Tx(:,Inod),neq,1);
        T=T+sparseUA(MUA.connectivity(:,Inod)+neqx,ones(MUA.Nele,1),Ty(:,Inod),neq,1);
        T=T+sparseUA(MUA.connectivity(:,Inod)+2*neqx,ones(MUA.Nele,1),Th(:,Inod),neq,1);
        
        F=F+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Fx(:,Inod),neq,1);
        F=F+sparseUA(MUA.connectivity(:,Inod)+neqx,ones(MUA.Nele,1),Fy(:,Inod),neq,1);
        F=F+sparseUA(MUA.connectivity(:,Inod)+2*neqx,ones(MUA.Nele,1),Fh(:,Inod),neq,1);
    end
    %%
    
    R=T-F;
    
    if ~Ronly
        
        
        
        % large memory version
        largeMemory=1;
        if largeMemory==1
            
            Iind=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1); Jind=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1);Xval=zeros(9*MUA.nod*MUA.nod*MUA.Nele,1);
            istak=0;
            for Inod=1:MUA.nod
                for Jnod=1:MUA.nod
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Kxu(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Kxv(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Kxh(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Kyu(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Kyv(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Kyh(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Khu(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Khv(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Khh(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                end
                
            end
            
            %%
            
            %            tSparse=tic;
            %            K=sparse(Iind,Jind,Xval,neq,neq);
            %            tSparse=toc(tSparse)   ;
            
            %            tSparse2=tic;
            K=sparseUA(Iind,Jind,Xval,neq,neq);
            %            tSparse2=toc(tSparse2)    ;
            
            
            %nz=nnz(K);
            
            
            %%
            
        else
            Iind=zeros(9*MUA.nod*MUA.Nele,1); Jind=zeros(9*MUA.nod*MUA.Nele,1);Xval=zeros(9*MUA.nod*MUA.Nele,1);
            K=sparseUA(neq,neq);
            
            for Inod=1:MUA.nod
                istak=0;
                for Jnod=1:MUA.nod
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Kxu(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Kxv(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Kxh(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Kyu(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Kyv(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Kyh(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=Khu(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=Khv(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                    
                    Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+MUA.Nele)=Khh(:,Inod,Jnod);
                    istak=istak+MUA.Nele;
                end
                K=K+sparseUA(Iind,Jind,Xval,neq,neq);
            end
        end
    end
    
    %    Kxu2=Kxu; Kxv2=Kxv; Kxh2=Kxh;  Kyu2=Kyu; Kyv2=Kyv;  Kyh2=Kyh ; Khu2=Khu; Khv2=Khv ; Khh2=Khh;
    %    save File2 Kxu2 Kxv2 Kxh2 Kyu2 Kyv2 Kyh2 Khu2 Khv2 Khh2
    
    if CtrlVar.IncludeTG3uvhBoundaryTerm && CtrlVar.TG3
        %[Ktest,Rtest]=BoundaryIntegralFullyImplicitTG3(CtrlVar,MUA,h0,h,u0,v0,u,v,as0+ab0,as1+ab1,dt);
        %[K,rh]=BoundaryIntegralFullyImplicitTG3(coordinates,connectivity,Boundary,h0,h1,u0,v0,u1,v1,a0,a1,dt,CtrlVar)
         [Ktest,Rtest]=BoundaryIntegralFullyImplicitTG3(MUA.coordinates,MUA.connectivity,MUA.Boundary,h0,h,u0,v0,u,v,as0+ab0,as1+ab1,dt,CtrlVar);
        R=R+Rtest;
        if ~Ronly
            K=K+Ktest;
        end
    end
    
    minh=min(h);
    

    if minh<2*CtrlVar.ThickMin && CtrlVar.InfoLevelNonLinIt>100   % if min thickness is approaching ThickMin give some information on h within NR loop
        msg=sprintf('In NRuvh loop, assembly stage: min(h) %-f \t max(h) %-g \n ',minh,max(h)) ;
        fprintf(CtrlVar.fidlog,msg) ; 
    end
    
    if ~Ronly  
        if full(any(isnan(diag(K))))
            save TestSave  ;  
            error(' NaN in K ' ) ;
        end ; 
    end
    
    if any(isnan(R))  
        save TestSave  ;  
        error(' NaN in R ' ) ; 
    end
end




