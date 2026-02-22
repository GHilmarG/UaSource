







function dIdA=dIdAq(CtrlVar,UserVar,MUA,F,BCS,BCsAdjoint,uAdjoint,vAdjoint,Meas)



%%
%
% Calculates
%
%
%
% $$ \langle  \delta_{A_i} F^x | \lambda_x \rangle + \langle  \delta_{A_i} F^y | \lambda_y \rangle $$
%
% where
%
%
% $$ F_x=\partial_x ( h \eta ( 4 \partial_x u + 2 \partial_y v)) + \partial_y ( h \eta (\partial_y u + \partial_x v) ) - t_x -   \frac{1}{2} g \partial_x (\rho h^2 -  \rho_o d^2)- g\,\mathcal{H}(h-h_f) (\rho h -\rho_o H^{+}) \partial_x B =0 $$
%
% $$ F_y = \partial_y ( h \eta ( 4 \partial_y v + 2 \partial_x u)) + \partial_x ( h \eta (\partial_x v + \partial_y y) ) - t_y -   \frac{1}{2} g \partial_y (\rho h^2 -  \rho_o d^2)- g\,\mathcal{H}(h-h_f) (\rho h -\rho_o H^{+}) \partial_y B =0 $$
% 
% with
%
% $$ t_x = t_x(C,u,v) $$
%
% $$ t_y = t_y(C,u,v) $$
%
% and
%
% $$ \eta=\eta(A,u,v) $$
%
% where 
%
% $$ \eta = \frac{1}{2} A^{-1/n} \, e^{(1-n)/n} +\eta_0 $$  
%
% with 
%
% $$ e=\sqrt{\epsilon_0^2+\epsilon_{xx}^2+\epsilon_{yy}^2+\epsilon_{xx} \, \epsilon_{yy}+\epsilon_{xy}^2} $$
%
% The Finite-Element terms of the momentum equations involving $\eta$ are
% 
% $$ 
% \langle  F_x | \phi_k \rangle =-\langle h \, \eta \, ( 4 \partial_x u + 2 \partial_y v)) | \partial_x \phi_k \rangle - \langle h \, \eta \, (\partial_y u + \partial_x v) | \partial_y \phi_k \rangle  
% $$
%
%
% $$ 
% \langle  F_y | \phi_k \rangle =-\langle h \, \eta \, ( 4 \partial_y v + 2 \partial_x u)) | \partial_y \phi_k \rangle - \langle h \, \eta \, (\partial_x v + \partial_y u) | \partial_x \phi_k \rangle  
% $$
%
% and 
%
%
% $$ \delta_A  \eta   = -\frac{1}{2n}  \, A^{-1/n-1} \; e^{(1-n)/n} \; \delta A $$
%
% or
%
% $$ \delta_A  \eta   = \partial_A \eta \; \delta A $$
%
%
% The second derivative is:
%
% $$ \delta^2_{AA}  \eta   = \frac{1}{2} \;  \frac{-1}{n} \; (-1/n-1)   \, A^{-1/n-2} \; e^{(1-n)/n} \; \delta A \, \delta A $$
%
% i.e.
%
% $$ \delta^2_{AA}  \eta   = \frac{1}{2} \;  (1/n^2 + 1/n)   \, A^{-1/n-2} \; e^{(1-n)/n} \; \delta A \, \delta A $$
%
% or
%
% $$ \delta^2_{AA}  \eta   = \partial^2_{AA} \eta  \; \delta A \, \delta A $$
%
% We therefore have
%
% $$ 
% \langle  \delta_{A_i} F^x | \lambda_x \rangle = -\langle h \, (d\eta/dA) \, \phi_i \, ( 4 \partial_x u + 2 \partial_y v) | \partial_x \lambda_x \rangle - \langle h \, (d \eta/dA) \, \phi_i \, (\partial_y u + \partial_x v) | \partial_y \lambda_x \rangle  
% $$
%
% that is
%
% $$ 
% \langle  \delta_{A_i} F^x | \lambda_x \rangle 
% = -\int \, \partial_A \eta \,  \, h \, \big (  ( 4 \partial_x u + 2 \partial_y v) \, \partial_x \lambda_x  + (\partial_y u + \partial_x v) \, \partial_y \lambda_x \big ) \, \phi_i \; dx \, dy 
% $$
%
% and
%
% $$ 
% \langle  \delta_{A_i} F^x | \lambda_x \rangle + \langle  \delta_{A_i} F^y | \lambda_y \rangle 
% = -\int \, \partial_A \eta \,  \, h \, \big (  ( 4 \partial_x u + 2 \partial_y v) \, \partial_x \lambda_x  + (\partial_y u + \partial_x v) \, \partial_y \lambda_x \big ) \, \phi_i \; dx \, dy 
%   -\int \, \partial_A \eta \,  \, h \, \big (  ( 4 \partial_y v + 2 \partial_x u) \, \partial_y \lambda_y  + (\partial_x v + \partial_y u) \, \partial_x \lambda_y \big ) \, \phi_i \; dx \, dy 
% $$
%
%
% which we can write as
%
% $$
% \langle  \delta_{A_i} F^x | \lambda_x \rangle + \langle  \delta_{A_i} F^y | \lambda_y \rangle 
% = -\int \, \partial_A \eta \,  \, h \, 
%   \big  (  ( 4 \partial_x u + 2 \partial_y v) \, \partial_x \lambda_x  + (\partial_y u + \partial_x v) \, \partial_y \lambda_x 
%  +         ( 4 \partial_y v + 2 \partial_x u) \, \partial_y \lambda_y  + (\partial_x v + \partial_y u) \, \partial_x \lambda_y \big ) \, \phi_i \; dx \, dy 
% $$
%
% This is a vector.
%
%
% $$
% \langle  \delta^2_{A_i\,A_j} F^x | \lambda_x \rangle + \langle  \delta^2_{A_i\,A_j} F^y | \lambda_y \rangle 
% = -\int \, \partial^2_{AA} \eta \,  \, h \, 
%   \big  (  ( 4 \partial_x u + 2 \partial_y v) \, \partial_x \lambda_x  + (\partial_y u + \partial_x v) \, \partial_y \lambda_x 
%  +         ( 4 \partial_y v + 2 \partial_x u) \, \partial_y \lambda_y  + (\partial_x v + \partial_y u) \, \partial_x \lambda_y \big ) \, \phi_i \, \phi_j \; dx \, dy 
% $$
%
% This is a matrix
%
%%








narginchk(9,9)

ndim=2;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
uAdjointnod=reshape(uAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
vAdjointnod=reshape(vAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
AGlennod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

% [points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);


for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;

    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        detJ=MUA.DetJ(:,Iint);
        Deriv=MUA.Deriv(:,:,:,Iint);
    else
        [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end

    hint=hnod*fun;
    nint=nnod*fun;
    AGlenInt=AGlennod*fun;
    AGlenInt(AGlenInt<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;


    exx=zeros(MUA.Nele,1); exy=zeros(MUA.Nele,1);  eyy=zeros(MUA.Nele,1);
    dlxdx=zeros(MUA.Nele,1); dlydx=zeros(MUA.Nele,1); dlxdy=zeros(MUA.Nele,1); dlydy=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod

        exx=exx+Deriv(:,1,Inod).*unod(:,Inod);
        eyy=eyy+Deriv(:,2,Inod).*vnod(:,Inod);
        exy=exy+0.5*(Deriv(:,1,Inod).*vnod(:,Inod) + Deriv(:,2,Inod).*unod(:,Inod));


        dlxdx=dlxdx+Deriv(:,1,Inod).*uAdjointnod(:,Inod);
        dlydx=dlydx+Deriv(:,1,Inod).*vAdjointnod(:,Inod);
        dlxdy=dlxdy+Deriv(:,2,Inod).*uAdjointnod(:,Inod);
        dlydy=dlydy+Deriv(:,2,Inod).*vAdjointnod(:,Inod);

    end


    detJw=detJ*MUA.weights(Iint);


    [~,~,~,dEtadA]=EffectiveViscositySSTREAM(CtrlVar,AGlenInt,nint,exx,eyy,exy);
    %dEtadA=dEtadA.*hint;




    for Inod=1:MUA.nod
        T(:,Inod)=T(:,Inod)...
            -dEtadA.*hint.*((4*exx+2*eyy).*dlxdx+2*exy.*dlxdy+(4*eyy+2*exx).*dlydy+2*exy.*dlydx).*fun(Inod).*detJw;
        %-dEtadA.*((4*exx+2*eyy).*dlxdx+(dudy+dvdx).*dlxdy+(4*eyy+2*exx).*dlydy+(dudy+dvdx).*dlydx).*fun(Inod).*detJw;



    end
end

dIdA=zeros(MUA.Nnodes,1);

for Inod=1:MUA.nod
    dIdA=dIdA+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end


if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
    dIdA=log(10)*F.AGlen.*dIdA;
end


dIdA=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCsAdjoint,dIdA);

end






