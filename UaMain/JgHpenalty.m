
function [J,g,H]=JgHpenalty(CtrlVar,MUA,x,x0,k,a)

%%
%
% *Terms involving functions of the latent variables*
%
% If we have a term involving a function of a variable, e.g.
%
% $$f(B)$$
%
% with
%
% $$
% R =\frac{1}{2} \int (f(B))^2 \, \mathrm{d}x
% = \frac{1}{2} \int (f(B)_i \phi_i(x) ) \, ( f(B)_j \phi_j(x) ) \, dx
% = \frac{1}{2} f_i(B) \left  (  \int \phi_i(x) ) \, \phi_j(x) \, dx \right ) \, f_j(B)
% = \frac{1}{2} f_i(B) \, M_{ij} \, f_j(B)
% $$
%
% where, for example, we could have
%
% $$ f(B)=\mathrm{SoftPlus}(B-s+h_{\mathrm{min}}) $$
%
%
% $$ f(B) = \sum_{k=1}^n \, f_k(B_k) \, \phi_k(x) $$
%
% Denote the derivative with $g$ , i.e.
%
% $$ g(B) := \partial f(B)/\partial B = \sum_{i=1}^n \, \frac{\partial f_k}{\partial B} \, \phi_k(x) = \sum_{k=1}^n \, g_k \, \phi_k(x) $$
%
% Note that in the case where, at every node
%
% $$g_k=\frac{\partial f_k}{\partial B}=1 $$
%
% for every $g_k$, we have
%
% $$
% g(B) := \partial f(B)/\partial B
% = \sum_{k=1}^n \, \frac{\partial f_k}{\partial B} \, \phi_k(x)
% = \sum_{k=1}^n \, g_k \, \phi_k(x)
% = \sum_{k=1}^n \, 1 \, \phi_k(x)  = 1
% $$
%
% due to the partition of unity property of the finite element form functions
%
% $$\sum_{k=1}^n \, \phi_k(x)  = 1 $$
%
%
% In general
%
% $$
% \delta_B R(B) = \lim_{\epsilon \to 0} \, \frac{1}{2} \frac{d}{d \epsilon} \int  (f(B(x)+\epsilon \phi_q(x)))^2 \, dx
% =\int f(B) \, g(B)\, \phi_q  \; dx
% = \int \; f_i \phi_i(x) \; g_k \phi_k(x) \; \phi_q(x)  \; dx
% $$
%
% In the above expression, we have taken the derivatives of $f(B)$ at the nodes and then interpolated those nodal values to the
% integration points. We could alternative interpolate $B$ to the integration points, and then take the derivative of $f(B)$
% with respect of $B$ at the location of the Integration points.
%
% Doing so results in:
%
% $$
% \delta_B R(B)
% = \int \; f_i \phi_i(x) \; \frac{\partial f(B_k \phi_k(x))}{\partial B}  \; \phi_q(x)  \; dx
% $$
%
% The term
%
% $$ \frac{\partial f(B_k \phi_k(x))}{\partial B}$$
%
% is a scalar, i.e. this derivative at a given location $x$
%
% Generally, when taking several derivatives and making use of the chain rule, it is best to evaluation derivative and other
% functions of the primary variables at the integration points.
%
% The Hessian is then
%
% $$
% H = \delta_{BB} I(B)
% = \int \; \frac{\partial f(B_i \phi_i(x)) }{\partial B} \, \phi_r(x) \; \frac{\partial f(B_k \phi_k(x))}{\partial B}  \; \phi_q(x)  \; dx
% + \int \; f_i \phi_i(x)  \; \frac{\partial^2 f(B_k \phi_k(x))}{\partial B \, \partial B}  \; \phi_r(x) \, \phi_q(x)  \; dx
% $$
%
% Note that for
%
% $$ f(B)= B = B_i \phi_i(x) $$
%
% in which case
%
% $$ f_i=B_I $$
%
% and
%
% $$f^{\prime} =1 $$
%
% and
%
% $$f^{\prime\,\prime} =0  $$
%
%
% we have
%
% $$ \frac{\partial p}{\partial B_j} = \delta_{ij} $$
%
%
% and we arrive at
%
% $$
% \delta_B R(B) = \lim_{\epsilon \to 0} \, \frac{1}{2} \frac{d}{d \epsilon} \int  (f(B(x)+\epsilon \phi_q(x)))^2 \, dx
% = f_i \; \left ( \int \phi_i(x) \, \phi_j(x) \, dx \right ) \frac{\partial f(B)}{\partial B_j}
% = f_i M_{ij} \phi_j \, \delta_{ij} = B_i M_{ij} = M_{ji} B_j
% $$
%
% or we then have (as listed above)
%
% $$ R= \frac{1}{2} \mathbf{B}^T \mathbf{M} \mathbf{B} $$
%
% $$\delta_B R(B) = \mathbf{M} \mathbf{B} $$
%
% and
%
% $$ H= \mathbf{M} $$
%
%
%%
%
%
%%

ndim=2;  neq=MUA.Nnodes;

xnod=reshape(x(MUA.connectivity,1),MUA.Nele,MUA.nod);
x0nod=reshape(x0(MUA.connectivity,1),MUA.Nele,MUA.nod);


%[points,weights]=sample('triangle',nip,ndim);

Jint=zeros(MUA.Nele,1); 
gint=zeros(MUA.Nele,MUA.Nnodes);
Hint=zeros(MUA.Nele,MUA.nod,MUA.nod);


for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    %Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
    detJw=detJ*MUA.weights(Iint);

    %% calculate J
    xint=xnod*fun;
    x0int=x0nod*fun;
 
    [f,dfdx,ddfdxdx] = SoftPlus(k,xint,x0int,Plot=false);

    % f=Bint; dfdB=1; ddfdBdB=0; % this results in ddfdBdB being equal to the mass matrix

    f2=a*f.*f/2;   % evaluated at integration point

    Jint=Jint+f2.*detJw;
    %%

   %% Calculate g


    for Inod=1:MUA.nod

        gint(:,Inod)=gint(:,Inod) + a*f.*dfdx.*fun(Inod).*detJw; % evaluated at integration point


        for Jnod=1:MUA.nod
            Hint(:,Inod,Jnod)=Hint(:,Inod,Jnod)...
                +a*(dfdx.*fun(Jnod).* dfdx.* fun(Inod)+f.*ddfdxdx.*fun(Jnod).*fun(Inod))...
            .*detJw;
        end

    end
end

J=sum(Jint);

ig=zeros(MUA.nod*MUA.Nele,1,"uint32");
One=ones(1,1,"uint32");
Xval=zeros(MUA.nod*MUA.Nele,1);
istak=0;

for Inod=1:MUA.nod


    ig(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
    Xval(istak+1:istak+MUA.Nele)=gint(:,Inod);
  
    istak=istak+MUA.Nele;

end

g=sparseUA(ig,One,Xval,neq,1);
g=full(g);

Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
istak=0;
for Inod=1:MUA.nod
    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Xval(istak+1:istak+MUA.Nele)=Hint(:,Inod,Jnod);
        istak=istak+MUA.Nele;

    end
end

H=sparseUA(Iind,Jind,Xval,neq,neq);
H=(H+H.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strictly so


 

end



