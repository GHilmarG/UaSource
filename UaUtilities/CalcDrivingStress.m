function [tx,ty]=CalcDrivingStress(CtrlVar,MUA,rho,g,s,h)


%  calculates (SIA) driving stress
%  tx=(rho g)^n | grad_{xy} s|^(n-1) h^(n+1) \p_x s
%  ty=(rho g)^n | grad_{xy} s|^(n-1) h^(n+1) \p_y s


hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);
snod=reshape(s(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

ndim=2;
points=sample('triangle',MUA.nip,ndim);
tx=zeros(MUA.Nele,MUA.nip);
ty=zeros(MUA.Nele,MUA.nip);

% vector over all elements for each integartion point
for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ')
        Deriv=MUA.Deriv(:,:,:,Iint);
    else
        Deriv=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    
    
    hint=hnod*fun;
    rhoint=rhonod*fun;
    
    
    dsdx=zeros(MUA.Nele,1); dsdy=zeros(MUA.Nele,1);
    
    
    for Inod=1:MUA.nod
        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
    end
    
    %gradSurf=sqrt(abs(dsdx.*dsdx+dsdy.*dsdy));
        
    T=rhoint.*g.*hint;
    tx(:,Iint)=-T.*dsdx;
    ty(:,Iint)=-T.*dsdy;
    
end

[tx,ty]=ProjectFintOntoNodes(MUA,tx,ty);


end