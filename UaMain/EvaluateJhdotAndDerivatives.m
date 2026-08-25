





function [Jhdot,duJhdot,dvJhdot,dhJhdot]=EvaluateJhdotAndDerivatives(UserVar,CtrlVar,MUA,F,BCs,Meas)

%%  Provides cost function and derivatives with respect to $u$, $v$, and $h$, of the cost function term involving $\dot{h}$
%
% $$J_{\dot{h}} = \| \dot{h} - \hat{\dot{h}} \| $$
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
% $$ J_{\dot{h}} = \frac{1}{2 \mathcal{A}} \int \! \int \left  ( \frac{\dot{h} - \hat{\dot{h}}}{h_{err}}  \right )^2 \; dx \, dy $$ 
%
% $$ \delta_u J_{\dot{h}} = \frac{1}{\mathcal{A}} \int \! \int \frac{\dot{h} - \hat{\dot{h}}}{h_{err}^2}  \, \delta_u \dot{h} \; dx \, dy $$ 
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
% $$ \delta_v \dot{h} = -\frac{1}{\rho} \partial_x ( \rho \, h \, \delta v )  $$
%
% $$ \delta_h \dot{h} = -\frac{1}{\rho} \partial_x ( \rho \, u \delta h + \rho \, v \, \delta h)  $$
%
% And therefore
%
% $$ \delta_u J_{\dot{h}} = -\frac{1}{\mathcal{A}} \int \! \int \frac{\dot{h} - \hat{\dot{h}}}{h_{err}^2}  \, \frac{1}{\rho} \partial_x ( \rho \, h \, \delta u )  \; dx \, dy $$ 
%
%
%
%
% $$
% r = \dot h_{obs} - a 
% + \left(\partial_x h+\frac{h}{\rho}\partial_x\rho\right)u + h\,\partial_x u 
% + \left(\partial_y h+\frac{h}{\rho}\partial_y\rho\right)v + h\,\partial_y v 
% $$
%
% $$
% \delta_uJ_{\dot h}[\phi_i] = 
% \frac{1}{\mathcal A}\int r\,\epsilon_{\dot h}^{-2}\left[\left(\partial_xh+\frac{h}{\rho}\partial_x\rho\right)\phi_i +
% h\,\partial_x\phi_i\right] dx\,dy 
% $$
%
%
% $$ \delta_vJ_{\dot h}[\phi_i] = 
% \frac{1}{\mathcal A}\int r\,\epsilon_{\dot h}^{-2}\left[\left(\partial_yh+\frac{h}{\rho}\partial_y\rho\right)\phi_i + h\,\partial_y\phi_i\right]dx\,dy$$
%
% see also: dhdtExplicit.m
%
%%


ndim=2; dof=1; neq=dof*MUA.Nnodes;

anod=reshape(F.as(MUA.connectivity,1),MUA.Nele,MUA.nod)+reshape(F.ab(MUA.connectivity,1),MUA.Nele,MUA.nod);
hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

[~,F.dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs) ; 

if ~isempty(F.dhdt) || ~isnan(F.dhdt)
    dhdtnod=reshape(F.dhdt(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    dhdtnod=nan;
end

dhdtMeasnod=reshape(Meas.dhdt(MUA.connectivity,1),MUA.Nele,MUA.nod);


JhdotIntSum=zeros(MUA.Nele,1);
duJhdotIntSum=zeros(MUA.Nele,MUA.nod);
dvJhdotIntSum=zeros(MUA.Nele,MUA.nod);
dhJhdotIntSum=zeros(MUA.Nele,MUA.nod);

Area=TriAreaTotalFE(MUA.coordinates,MUA.connectivity);
dhdtErr=sqrt(spdiags(Meas.dhdtCov));
dhdtErrnod=reshape(dhdtErr(MUA.connectivity,1),MUA.Nele,MUA.nod);


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
    hdotErrInt=dhdtErrnod*fun;

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


     R=(hdot-hdotMeasint)./hdotErrInt; 


    detJw=detJ*MUA.weights(Iint);
    
  
    JhdotIntSum=JhdotIntSum+((hdot-hdotMeasint)./hdotErrInt).^2 .*detJw/2/Area; 
    
    for Inod=1:MUA.nod
        
        % hdot=aint-(dhdx.*uint+hint.*dudx +dhdy.*vint+hint.*dvdy) ; 
        
       

        duJhdotInt=-R...
            .*(drhodx.*hint.*fun(Inod)+rhoint.*dhdx.*fun(Inod)+rhoint.*hint.*Deriv(:,1,Inod))./(hdotErrInt.*rhoint)...
            .*detJw/Area;
        
        
        dvJhdotInt=-R...
            .*(drhody.*hint.*fun(Inod)+rhoint.*dhdy.*fun(Inod)+rhoint.*hint.*Deriv(:,2,Inod))./(hdotErrInt.*rhoint)...
            .*detJw/Area;
        
        dhJhdotInt=-R...
            .*( ...
             drhodx.*uint.*fun(Inod)+rhoint.*dudx.*fun(Inod)+rhoint.*uint.*Deriv(:,1,Inod)...
            +drhody.*vint.*fun(Inod)+rhoint.*dvdy.*fun(Inod)+rhoint.*vint.*Deriv(:,2,Inod) ...
            )./(hdotErrInt.*rhoint)...
            .*detJw/Area;

        
        
        duJhdotIntSum(:,Inod)=duJhdotIntSum(:,Inod)+duJhdotInt;
        dvJhdotIntSum(:,Inod)=dvJhdotIntSum(:,Inod)+dvJhdotInt;
        dhJhdotIntSum(:,Inod)=dhJhdotIntSum(:,Inod)+dhJhdotInt;
        
        
    end
end

% assemble right-hand side

Jhdot=sum(JhdotIntSum) ;

duJhdot=sparseUA(neq,1);
dvJhdot=sparseUA(neq,1);
dhJhdot=sparseUA(neq,1);
for Inod=1:MUA.nod
    
    duJhdot=duJhdot+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),duJhdotIntSum(:,Inod),neq,1);
    dvJhdot=dvJhdot+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),dvJhdotIntSum(:,Inod),neq,1);
    dhJhdot=dhJhdot+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),dhJhdotIntSum(:,Inod),neq,1);
    
end


duJhdot=full(duJhdot);
dvJhdot=full(dvJhdot);
dhJhdot=full(dhJhdot);


%% If F.dhdt is available this should give the same answer
% dhdtErr=sqrt(spdiags(Meas.dhdtCov)) ;  dhdtres=(F.dhdt-Meas.dhdt)./dhdtErr ;  JhdotTest=full(dhdtres'*MUA.M*dhdtres)/2/Area;
%%


% Don't apply this here!  Because duJhdot and dvJhdot contribute to the right-hand side of the Adjoint equations.
% However, dhJhdot does not and this derivative needs to be projected, but do this later
%
% [duJhdot,dvJhdot,dhJhdot]=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCs,duJhdot,dvJhdot,dhJhdot);



end