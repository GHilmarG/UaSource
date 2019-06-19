function dIdAGlendata=dIdAEleSteps(CtrlVar,MUA,lx,ly,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF)

% calculates dIdAGlendata for step shifts within each element, returning element values
%
     
    
    ndim=2;
    
    
    hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    unod=reshape(ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
    vnod=reshape(vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
    lxnod=reshape(lx(MUA.connectivity,1),MUA.Nele,MUA.nod);
    lynod=reshape(ly(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    %if ~CtrlVar.AGlenisElementBased
    %    AGlennod=reshape(AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
    %end
    
    [points,weights]=sample('triangle',MUA.nip,ndim);
    [~,~,~,~,~,~,~,e]=calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n);
    
    dIdAGlendata=zeros(MUA.Nele,1); EleArea=zeros(MUA.Nele,1);
    
    for Iint=1:MUA.nip
        
        
        fun=shape_fun(Iint,ndim,MUA.nod,points) ;
        
        if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
            Deriv=MUA.Deriv(:,:,:,Iint);
            detJ=MUA.DetJ(:,Iint);
        else
            [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
        end
        
        hint=hnod*fun;
        %if ~CtrlVar.AGlenisElementBased
        %    AGlenInt=AGlennod*fun;
        %else
            AGlenInt=AGlen;
        %end
        
        dudx=zeros(MUA.Nele,1); dvdx=zeros(MUA.Nele,1); dudy=zeros(MUA.Nele,1); dvdy=zeros(MUA.Nele,1);
        dlxdx=zeros(MUA.Nele,1); dlydx=zeros(MUA.Nele,1); dlxdy=zeros(MUA.Nele,1); dlydy=zeros(MUA.Nele,1);
        
        for Inod=1:MUA.nod  % derivaties at integration points
            
            dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
            dvdx=dvdx+Deriv(:,1,Inod).*vnod(:,Inod);
            dudy=dudy+Deriv(:,2,Inod).*unod(:,Inod);
            dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
            
            dlxdx=dlxdx+Deriv(:,1,Inod).*lxnod(:,Inod);
            dlydx=dlydx+Deriv(:,1,Inod).*lynod(:,Inod);
            dlxdy=dlxdy+Deriv(:,2,Inod).*lxnod(:,Inod);
            dlydy=dlydy+Deriv(:,2,Inod).*lynod(:,Inod);
            
        end
        
        detJw=detJ*weights(Iint);
        EleArea=EleArea+detJw;
        
        
        dEtadA=-real(hint.*(AGlenInt+CtrlVar.AGlenAdjointZero).^(-1./n-1).*(e(:,Iint)+CtrlVar.AdjointEpsZero).^((1-n)./n))./(2*n);
        %dEtadA=-real(hint.* (AGlenInt.*e(:,Iint)+CtrlVar.AGlenAdjointZero).^(-(1+n)/n)  .*e(:,Iint).^(2/n)) /(2*n);
        
        if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
            
            dEtadA=log(10)*AGlenInt.*dEtadA;
            
        end
        dIdAGlendata=dIdAGlendata-dEtadA.*((4*dudx+2*dvdy).*dlxdx+(dudy+dvdx).*dlxdy+(4*dvdy+2*dudx).*dlydy+(dudy+dvdx).*dlydx).*detJw;
    end
    
    
    
    switch CtrlVar.Inverse.AdjointGradientPreMultiplier
        
        case 'M'
            
            
            if CtrlVar.Inverse.InfoLevel>=1000
                figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdAGlendata) ; title('dIdA Mesh Dependend')
            end
            % make mesh independent by dividing with element areas
            dIdAGlendataNorm=norm(dIdAGlendata);
            dIdAGlendata=dIdAGlendata./EleArea;
            dIdAGlendata=dIdAGlendata*dIdAGlendataNorm/norm(dIdAGlendata);
            
            
            if CtrlVar.Inverse.InfoLevel>=1000
                figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdAGlendata) ; title('dIdA Mesh Independend')
            end
            
    end
    
    
end



