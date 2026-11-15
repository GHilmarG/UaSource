





function [Jhdot,duJhdot,dvJhdot,dhJhdot]=EvaluateJhdotAndDerivatives(UserVar,CtrlVar,MUA,F,BCs,Meas)

%%  Provides cost function and derivatives with respect to $u$, $v$, and $h$, of the cost function term involving $\dot{h}$
%
% 
%
% $$ J_{\dot{h}} = \frac{1}{2 \mathcal{A}} \int \! \int \left  ( \frac{\dot{h} - \hat{\dot{h}}}{\dot{h}_{err}}  \right )^2 \; dx \, dy $$ 
%
%
% where
%
% $$ \rho \, \dot{h}  = \rho \, a -  \nabla \cdot \mathbf{q} $$
%
% or
%
% $$ \dot{h}  = a -  \frac{1}{\rho} \nabla \cdot \mathbf{q} $$
%
% Here
%
% $$\nabla \cdot \mathbf{q} = \partial_x (\rho h u ) + \partial_y (\rho h v) $$
%
% and
%
% $$ J_{\dot{h}} = \frac{1}{2 \mathcal{A}} \int \! \int \left  ( \frac{\dot{h} - \hat{\dot{h}}}{\dot{h}_{err}}  \right )^2 \; dx \, dy $$ 
%
% $$ \delta_u J_{\dot{h}} = \frac{1}{\mathcal{A}} \int \! \int \frac{\dot{h} - \hat{\dot{h}}}{\dot{h}_{err}^2}  \, \delta_u \dot{h} \; dx \, dy $$ 
%
% For
%
% $$ \dot{h}  = a -  \frac{1}{\rho} \nabla \cdot \mathbf{q} $$
%
% we find
%
% $$ 
% \delta_u \dot{h} = \lim_{\epsilon \to 0} \frac{d}{d\epsilon} ( a- \frac{1}{\rho} (\partial_x (\rho \, h \, ( u+ \epsilon \delta u ) + \partial_y ( \rho \, h \, v ) ) 
% = -\frac{1}{\rho} \partial_x ( \rho \, h \, \delta u ) 
% $$
% 
% that is
%
% $$ \delta_u \dot{h} = -\frac{1}{\rho} \partial_x ( \rho \, h \, \delta u )  $$
%
% $$ \delta_v \dot{h} = -\frac{1}{\rho} \partial_y ( \rho \, h \, \delta v )  $$
%
% $$ \delta_h \dot{h} = -\frac{1}{\rho} (\partial_x ( \rho \, u \delta h ) + \partial_y ( \rho \, v \, \delta h) ) $$
%
% And therefore
%
% $$ \delta_u J_{\dot{h}} = -\frac{1}{\mathcal{A}} \int \! \int \frac{\dot{h} - \hat{\dot{h}}}{\dot{h}_{err}^2}  \, \frac{1}{\rho} \partial_x ( \rho \, h \, \delta u )  \; dx \, dy $$ 
%
%
%
% $$ R = \frac{\dot{h} - \dot{h}_{obs}}{\dot{h}_{err}} $$ 
%
% $$
% \delta_uJ_{\dot h}[\phi_i] = 
% -\frac{1}{\mathcal A}\int R\,\dot{h}_{err}^{-1}\left[\left(\partial_xh+\frac{h}{\rho}\partial_x\rho\right)\phi_i +
% h\,\partial_x\phi_i\right] dx\,dy 
% $$
%
%
% $$ \delta_vJ_{\dot h}[\phi_i] = 
% -\frac{1}{\mathcal A}\int R\,\dot{h}_{err}^{-1}\left[\left(\partial_yh+\frac{h}{\rho}\partial_y\rho\right)\phi_i + h\,\partial_y\phi_i\right]dx\,dy$$
%
% see also: dhdtExplicit.m
%
% Assumptions: It is here assumed that $a=a_b+a_s$ are independent of $u$, $v$, and $h$. For if, for example, the basal
% melt-rate were to be described as a function of draft, that case would not be correctly treated..
%
%
%
%%


%% "-points-" and "-blocks-" formulations (CtrlVar.Inverse.DataMisfit.dhdt.Method), see dhdtObservationMisfit.m (14 November 2026)
%
% In these formulations the modelled dh/dt is brought to the data: evaluated at the measurement locations, or averaged over the blocks
% of the measurements. Errors are treated as uncorrelated, and there is no normalisation with the area of the domain.
if dhdtMisfitMethod(CtrlVar)~="-nodal-integral-"
    [Jhdot,duJhdot,dvJhdot,dhJhdot]=dhdtObservationMisfit(CtrlVar,MUA,F,Meas);
    return
end

%% "-nodal-integral-" formulation (default)

ndim=2; dof=1; neq=dof*MUA.Nnodes;

anod=reshape(F.as(MUA.connectivity,1),MUA.Nele,MUA.nod)+reshape(F.ab(MUA.connectivity,1),MUA.Nele,MUA.nod);
hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

%[~,F.dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs) ;



dhdtMeasnod=reshape(Meas.dhdt(MUA.connectivity,1),MUA.Nele,MUA.nod);


JhdotIntSum=zeros(MUA.Nele,1);
duJhdotIntSum=zeros(MUA.Nele,MUA.nod);
dvJhdotIntSum=zeros(MUA.Nele,MUA.nod);
dhJhdotIntSum=zeros(MUA.Nele,MUA.nod);

Area=TriAreaTotalFE(MUA.coordinates,MUA.connectivity);
% Weights w=1/dhdtErr^2 are defined at the nodes and interpolated to the integration points. This is the same discretisation
% as used for the velocity misfit term in Misfit.m (weighted mass matrix with nodal weights w), so that both misfit terms are
% treated consistently. (Previously the error itself was interpolated and then squared.)
%
% Note: diag() rather than spdiags(). With one output, spdiags returns all non-zero diagonals as columns, which would give
% wrong results should dhdtCov ever contain off-diagonal entries. full() avoids sparse arithmetic in the integration loop.
dhdtW=1./full(diag(Meas.dhdtCov));

dhdtWnod=reshape(dhdtW(MUA.connectivity,1),MUA.Nele,MUA.nod);


% vector over all elements for each integration point
for Iint=1:MUA.nip


    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);

    aint=anod*fun;
    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    rhoint=rhonod*fun;


    hdotMeasint=dhdtMeasnod*fun;
    hdotWInt=dhdtWnod*fun;    % weight, 1/dhdtErr^2, at integration points

    dhdx=zeros(MUA.Nele,1);
    dhdy=zeros(MUA.Nele,1);
    dudx=zeros(MUA.Nele,1);
    dvdy=zeros(MUA.Nele,1);
    drhodx=zeros(MUA.Nele,1);
    drhody=zeros(MUA.Nele,1);
    % derivatives at one integration point for all elements
    for Inod=1:MUA.nod

        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);

        dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
        dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);

        drhodx=drhodx+Deriv(:,1,Inod).*rhonod(:,Inod);
        drhody=drhody+Deriv(:,2,Inod).*rhonod(:,Inod);
    end


   % for this to be a consistent derivative, I must evaluate this directly from the integration point values.
   % Although F.dhdt has already been calculated by projecting onto the nodes, I must here use hdot calculated at int 
   %
   % if ~isnan(dhdtnod)
   %     hdot=dhdtnod*fun;
   % else
        hdot=aint-(rhoint.*dhdx.*uint+rhoint.*hint.*dudx+drhodx.*hint.*uint+rhoint.*dhdy.*vint+rhoint.*hint.*dvdy+drhody.*hint.*vint)./rhoint ;
    % end


     r=hdot-hdotMeasint;   % residual at integration points


    detJw=detJ*MUA.weights(Iint);
    
  
    JhdotIntSum=JhdotIntSum+hdotWInt.*r.^2 .*detJw/2/Area;
    
    for Inod=1:MUA.nod
        
        % hdot=aint-(dhdx.*uint+hint.*dudx +dhdy.*vint+hint.*dvdy) ; 
        
       

        duJhdotInt=-hdotWInt.*r...
            .*(drhodx.*hint.*fun(Inod)+rhoint.*dhdx.*fun(Inod)+rhoint.*hint.*Deriv(:,1,Inod))./rhoint...
            .*detJw/Area;
        
        
        dvJhdotInt=-hdotWInt.*r...
            .*(drhody.*hint.*fun(Inod)+rhoint.*dhdy.*fun(Inod)+rhoint.*hint.*Deriv(:,2,Inod))./rhoint...
            .*detJw/Area;
        
        dhJhdotInt=-hdotWInt.*r...
            .*( ...
             drhodx.*uint.*fun(Inod)+rhoint.*dudx.*fun(Inod)+rhoint.*uint.*Deriv(:,1,Inod)...
            +drhody.*vint.*fun(Inod)+rhoint.*dvdy.*fun(Inod)+rhoint.*vint.*Deriv(:,2,Inod) ...
            )./rhoint...
            .*detJw/Area;

        
        
        duJhdotIntSum(:,Inod)=duJhdotIntSum(:,Inod)+duJhdotInt;
        dvJhdotIntSum(:,Inod)=dvJhdotIntSum(:,Inod)+dvJhdotInt;
        dhJhdotIntSum(:,Inod)=dhJhdotIntSum(:,Inod)+dhJhdotInt;
        
        
    end
end

% assemble right-hand side

Jhdot=full(sum(JhdotIntSum)) ;

duJhdot=sparse(neq,1);
dvJhdot=sparse(neq,1);
dhJhdot=sparse(neq,1);
for Inod=1:MUA.nod
    
    duJhdot=duJhdot+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),duJhdotIntSum(:,Inod),neq,1);
    dvJhdot=dvJhdot+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),dvJhdotIntSum(:,Inod),neq,1);
    dhJhdot=dhJhdot+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),dhJhdotIntSum(:,Inod),neq,1);
    
end


duJhdot=full(duJhdot);
dvJhdot=full(dvJhdot);
dhJhdot=full(dhJhdot);


%% If F.dhdt is available this should give (approx) the same answer
% dhdtW=1./full(diag(Meas.dhdtCov)) ; [~,Mw]=GeneralizedMassAndStiffnessAssembly([],CtrlVar,MUA,dhdtW,0,0,0) ;
% dhdtres=F.dhdt-Meas.dhdt ;  JhdotTest=full(dhdtres'*Mw*dhdtres)/2/Area;
%%



end