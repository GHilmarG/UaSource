


function HPsiddFdCdC=PsiTimesddFuvdCdC(CtrlVar,MUA,F,uAdjoint,vAdjoint)

narginchk(5,5)

%%
%
% Calculates
%
% $$\Psi_n \frac{\partial^2 F_n}{\partial C_i \, \partial C_j} $$
%
%
% $$ \langle \Psi(x) \mid \delta^2_{CC} F(x) \rangle = \langle \Psi(x) \mid \partial^2_{CC} F(x) \, \delta C \, \delta C \rangle $$
%
% Or in a discrete form:
%
% $$\Psi_n \frac{\partial^2 F_n}{\partial C_i \, \partial C_j} $$
%
% For example, for Weertman
%
% $$F_x= -t_x=-\mathcal{G} \beta^2 \, v_x $$
%
% $$F_y=- t_y=-\mathcal{G} \beta^2 \, v_y $$
%
%
% $$ \beta^2 = (C+C_0)^{-1/m} \, (v_x^2+v_y^2 + v_0)^{(1-m)/2m} $$
%
% $$F_x= -\mathcal{G} \; (C+C_0)^{-1/m} \;  (v_x^2+v_y^2+v_0^2)^{(1-m)/2m}  \, v_x $$
%
% we have
%
% $$\delta_C F_x = \frac{\mathcal{G}}{m} \,  (C+C_0)^{-1/m-1} \; v^{1/m-1} \, v_x \; \delta C$$
%
% and
%
% $$\delta^2_{CC} F_x = \mathcal{G} \, \frac{1}{m} (-1/m-1)  \;  (C+C_0)^{-1/m-1} \; v^{1/m-1} \, v_x \; \delta C$$
%
% or 
%
% $$\delta^2_{CC} F_x = -\mathcal{G}  (1/m^2+1/m)  \;  (C+C_0)^{-1/m-1} \; v^{1/m-1} \, v_x \; \delta C$$
%
%
% $$
% \langle \lambda |  \delta_C F \rangle = \int \mathcal{G} \; \frac{1}{m}  \;  (C+C_0)^{-1/m-1} \; \delta C \; v^{1/m-1} \, \left (  \lambda_x \, v_x + \lambda_y \, v_y \right ) \; dx \, dy
% $$
%
% Which is a vector, and
%
% $$ \langle \Psi | \delta^2_{CC} F \rangle =  
%  -\int  \mathcal{G} \, (1/m^2+1/m) \;   (C+C_0)^{-1/m-2} \;  v^{1/m-1} \, \left (  \lambda_x \, v_x +
% \lambda_y \, v_y \right ) \; \delta C_i \, \delta C_j \; dx \, dy $$
% 
% which is a matrix.
%
% Note: I tested doing the C to log(C) change of variables either within the int loop, or afterwards on nodes. As far as I
% could see, the reduction in the Newton step was identical.
%
%
%%

if CtrlVar.SlidingLaw~="Weertman" 

    error("PsiTimesddFuvdCdC:NotImplemented","only implemented for Weertman sliding law.")

end

ndim=2; 
nNodes=MUA.Nnodes ;

C0=CtrlVar.Cmin;
u0=CtrlVar.SpeedZero;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod

Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
hfnod=F.rhow*(Snod-Bnod)./rhonod;

unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

uAdjointnod=reshape(uAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
vAdjointnod=reshape(vAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);

Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);



H=zeros(MUA.Nele,MUA.nod,MUA.nod);


for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    % Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);

    h=hnod*fun;   
    hf=hfnod*fun;
    
    u=unod*fun;
    v=vnod*fun;

    U=sqrt(u.*u + v.*v+u0^2) ; 

    C=Cnod*fun; 
    C(C<C0)=C0;

    m=mnod*fun;
    
    lx=uAdjointnod*fun;
    ly=vAdjointnod*fun;

    G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);

    detJw=detJ*MUA.weights(Iint);
   
    %Temp=  (1./(m+1)) .*  G.*  ((C+C0).^(-1./m-2)) .* ( U.^(1./m-1)) .*(lx.*u+ly.*v);  % at nodes

    Temp= (1./(m.^2) + 1./m) .*  G.*  ((C+C0).^(-1./m-2)) .* ( U.^(1./m-1)) .*(lx.*u+ly.*v).* (log(10).*C).^2; % at int

    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod

            H(:,Inod,Jnod)=H(:,Inod,Jnod) - Temp .*fun(Inod) .*fun(Jnod).*detJw;  % not sure I understand the sign here, 
            
        end
    end
end

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32'); 
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

istak=0;

for Inod=1:MUA.nod
    %istak=0;

    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); 
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); 
        Xval(istak+1:istak+MUA.Nele)=H(:,Inod,Jnod);
        %Xval(istak+1:istak+MUA.Nele)=H(:,Inod,Jnod).* (log(10).*F.C(MUA.connectivity(:,Inod))) .*(log(10).*F.C(MUA.connectivity(:,Jnod))) ; % At nodes
        
        istak=istak+MUA.Nele;

    end
end


HPsiddFdCdC=sparseUA(Iind,Jind,Xval,nNodes,nNodes);

HPsiddFdCdC=(HPsiddFdCdC'+HPsiddFdCdC)/2 ; 




end







