
function [exx,eyy,exy,e]=CalcHorizontalNodalStrainRates(CtrlVar,MUA,u,v)


%%
%  
% Calculates horizontal strain rates given horizontal velocities. 
%
%   [exx,eyy,exy,e]=CalcNodalStrainRates(CtrlVar,MUA,u,v)
%
% Returns nodal values.
%
% Note: The strain rates are calculated at integration points and then projected
% onto nodes. 
%
% The projection does not conserve positivity and positve integration
% values can become negative at nodes. The effectiv strain rate , e, is for
% this reason calculated directly from nodal values, ensuring that e is
% always positive.
%
%

unod=reshape(u(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(v(MUA.connectivity,1),MUA.Nele,MUA.nod);
exx=zeros(MUA.Nele,MUA.nip);
eyy=zeros(MUA.Nele,MUA.nip);
exy=zeros(MUA.Nele,MUA.nip);

for Iint=1:MUA.nip
    

    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        Deriv=MUA.Deriv(:,:,:,Iint);
       % detJ=MUA.DetJ(:,Iint);
    else
        Deriv=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,MUA.points,Iint);
    end

    
    for I=1:MUA.nod
        exx(:,Iint)=exx(:,Iint)+Deriv(:,1,I).*unod(:,I);
        eyy(:,Iint)=eyy(:,Iint)+Deriv(:,2,I).*vnod(:,I);
        exy(:,Iint)=exy(:,Iint)+0.5*(Deriv(:,1,I).*vnod(:,I) + Deriv(:,2,I).*unod(:,I));
    end
end


[exx,eyy,exy]=ProjectFintOntoNodes(MUA,exx,eyy,exy);
 
e=real(sqrt(exx.^2+eyy.^2+exx.*eyy+exy.^2));

end
