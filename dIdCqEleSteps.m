function dIdCdata=dIdCqEleSteps(CtrlVar,MUA,lx,ly,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF)
    % [dIdCdata]=dIdCqEleSteps(CtrlVar,MUA,ub,vb,ud,vd,lx,ly,S,B,h,C,m,rho,rhow,GF);
    %   
    %
    % calculates dIdC as an elementwise integration, returning element values
    %
        
    %
    %  C= C_q F_q ; where F_q(x,y)=1 for   (x,y) within element q, 0 otherwise
    %
    
    
    
    ndim=2;
    
    
    hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
    unod=reshape(ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
    vnod=reshape(vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
    Bnod=reshape(B(MUA.connectivity,1),MUA.Nele,MUA.nod);
    Snod=reshape(S(MUA.connectivity,1),MUA.Nele,MUA.nod);
    rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
    lxnod=reshape(lx(MUA.connectivity,1),MUA.Nele,MUA.nod);
    lynod=reshape(ly(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    
    if ~CtrlVar.CisElementBased
        Cnod=reshape(C(MUA.connectivity,1),MUA.Nele,MUA.nod);
    end
    
    
    [points,weights]=sample('triangle',MUA.nip,ndim);
    
    
    dIdCdata=zeros(MUA.Nele,1); EleArea=zeros(MUA.Nele,1);
    
    for Iint=1:MUA.nip
        
        
        fun=shape_fun(Iint,ndim,MUA.nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        %[~,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
        
         if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
            %Deriv=MUA.Deriv(:,:,:,Iint);
            detJ=MUA.DetJ(:,Iint);
        else
            [~,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
        end
        
        
        hint=hnod*fun;
        uint=unod*fun;
        vint=vnod*fun;
        
        if CtrlVar.CisElementBased==1
            Cint=C;
        else
            Cint=Cnod*fun;
        end

        %Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin; % for higher order elements it is possible that Cint is less than any of the nodal values
        Bint=Bnod*fun;
        Sint=Snod*fun;
        rhoint=rhonod*fun;
        lxint=lxnod*fun;
        lyint=lynod*fun;
        
        hfint=(Sint-Bint)*rhow./rhoint;
        Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0); 
        
        
        %Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
        Ctemp= -(1/m)*Heint.*(Cint+CtrlVar.CAdjointZero).^(-1/m-1).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1/m-1) ; % this is correct
        %Ctemp= -(1/m)*(Cint+CtrlVar.CAdjointZero).^(-1/m-1).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1/m-1) ;  % used for testing purposes
        
        
        detJw=detJ*weights(Iint);
        EleArea=EleArea+detJw;
        
        dIdCdata=dIdCdata-Ctemp.*(uint.*lxint+vint.*lyint).*detJw;
        
    end
    
    
    if CtrlVar.NormalizeWithAreas
        %dIdCNorm=norm(dIdC);
        dIdCdata=dIdCdata./EleArea;
        %dIdC=sum(EleArea)*dIdC;
        %dIdC=dIdC*norm(dIdCNorm)/norm(dIdC) ; % ensure that both estimates of dIdC have the same norm
    end
    
    
    
end



