function dIdC=dIdCqEleSteps(CtrlVar,MUA,uAdjoint,vAdjoint,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,GF)
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
uAdjointnod=reshape(uAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
vAdjointnod=reshape(vAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);


if ~CtrlVar.CisElementBased
    Cnod=reshape(C(MUA.connectivity,1),MUA.Nele,MUA.nod);
    mnod=reshape(m(MUA.connectivity,1),MUA.Nele,MUA.nod);
end


[points,weights]=sample('triangle',MUA.nip,ndim);
dIdC=zeros(MUA.Nele,1); EleArea=zeros(MUA.Nele,1);

for Iint=1:MUA.nip
    
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ;
    detJ=MUA.DetJ(:,Iint);
    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    
    if CtrlVar.CisElementBased==1
        Cint=C;
        mint=m;
    else
        Cint=Cnod*fun;
        Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
        mint=mnod*fun;
    end
    
    Bint=Bnod*fun;
    Sint=Snod*fun;
    rhoint=rhonod*fun;
    uAdjointint=uAdjointnod*fun;
    vAdjointint=vAdjointnod*fun;
    hfint=(Sint-Bint)*rhow./rhoint;
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    
    Ctemp= (1./mint).*Heint.*(Cint+CtrlVar.CAdjointZero).^(-1./mint-1).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1./mint-1) ;
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        Ctemp=log(10)*Cint.*Ctemp;
    end
    
    
    
    detJw=detJ*weights(Iint);
    EleArea=EleArea+detJw;
    
    dIdC=dIdC+Ctemp.*(uint.*uAdjointint+vint.*vAdjointint).*detJw;
    
end

switch  CtrlVar.Inverse.AdjointGradientPreMultiplier
    
    case 'M'
            
        if CtrlVar.Inverse.InfoLevel>=1000
            figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdC) ; title('dIdC Mesh Dependend')
        end
        
        % make mesh independent by dividing with element areas
        dIdCnorm=norm(dIdC);
        dIdC=dIdC./EleArea;
        dIdC=dIdCnorm*dIdC/norm(dIdC);
        
        if CtrlVar.Inverse.InfoLevel>=10
            fprintf('Making dIdC mesh independent by dividing with element areas.\n')
        end
             
        if CtrlVar.Inverse.InfoLevel>=1000
            figure ; PlotMeshScalarVariable(CtrlVar,MUA,dIdC) ; title('dIdC Mesh Independend')
        end
        
        
end


end



