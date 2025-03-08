




function R=ElementErrorEstimator(CtrlVar,MUA,F)

%
%
% Calculate field gradients at integration points (IP), and integrate over
% elements.
%
% Then L2 project the IP gradients onto nodal points (NP), and the integrate those over elements.
%


eInt=StrainRatesInt(CtrlVar,MUA,F);
eNode=ProjectFintOntoNodes(MUA,eInt);

eIntEle=FEintegrate2D(CtrlVar,MUA,eInt);
eNodeEle=FEintegrate2D(CtrlVar,MUA,eNode);


R=abs(eIntEle-eNodeEle);

end




function [e,exx,eyy,exy]=StrainRatesInt(CtrlVar,MUA,F)

ubnod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vbnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
coox=reshape(MUA.coordinates(MUA.connectivity,1),MUA.Nele,MUA.nod);
cooy=reshape(MUA.coordinates(MUA.connectivity,2),MUA.Nele,MUA.nod);

exx=zeros(MUA.Nele,MUA.nip); eyy=zeros(MUA.Nele,MUA.nip); exy=zeros(MUA.Nele,MUA.nip);
xint=zeros(MUA.Nele,MUA.nip) ; yint=zeros(MUA.Nele,MUA.nip);

ndim=2; 
for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;

    xint(:,Iint)=coox*fun;
    yint(:,Iint)=cooy*fun;

    Deriv=MUA.Deriv(:,:,:,Iint);

    for I=1:MUA.nod
        exx(:,Iint)=exx(:,Iint)+Deriv(:,1,I).*ubnod(:,I);
        eyy(:,Iint)=eyy(:,Iint)+Deriv(:,2,I).*vbnod(:,I);
        exy(:,Iint)=exy(:,Iint)+0.5*(Deriv(:,1,I).*vbnod(:,I) + Deriv(:,2,I).*ubnod(:,I));
    end
end

e=real(sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2));

end