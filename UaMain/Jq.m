



function [dJduv,Iuv,Ihdot]=Jq(CtrlVar,MUA,F,BCs,Meas)

narginchk(5,5)


%%
%
% Calculates the $$J^{q}$$ term, where $q$ is $u$ and $v$.
%
% There are number of terms in the cost function $J$, but the only ones which are explicit functions of $u$ and $v$ are:
%
% The $u$ misfit term:
%
% $$J_{u} = \frac{1}{2}\left(u_{\mathrm{obs}}-f_{u}\right)^{T}C_u^{-1}\left(u_{\mathrm{obs}}-f_v\right) $$
%
% The $v$ misfit term:
%
% $$J_{v} = \frac{1}{2}\left(v_{\mathrm{obs}}-f_{v}\right)^{T}C_v^{-1}\left(v_{\mathrm{obs}}-f_v\right) $$
%
% and the $\dot{h}$ misfit term:
%
% $$J_{\dot{h}}= \frac{1}{2\mathcal{A}} \int \left ( \dot{h}_{\mathrm{obs}} - a + \frac{1}{\rho} \left ( \partial_x (\rho h u)
% + \partial_y (\rho h v) \right ) \right )^2  \epsilon_{\dot{h}}^{-2}(x,y) \; dx \, dy
% $$
%
% Derivation and some useful quantities:
%
%
% $$
% r = \dot h_\mathrm{obs} - a + \left(\partial_xh+\frac{h}{\rho}\partial_x\rho\right)u + h\,\partial_xu + \left(\partial_yh+\frac{h}{\rho}\partial_y\rho\right)v + h\,\partial_yv
% $$
%
%
% $$
% \delta^2_{uu}J_{\dot h}[\phi_i,\phi_k] = \frac{1}{\mathcal A}\int\Big(\delta_ur[\phi_i]\Big)\Big(\delta_ur[\phi_k]\Big)\epsilon_{\dot h}^{-2}\,dx\,dy
% + \underbrace{\frac{1}{\mathcal A}\int r\,\delta^2_{uu}r[\phi_i,\phi_k]\,\epsilon_{\dot h}^{-2}\,dx\,dy}_{=0,\ \mathrm{since}\;\delta^2_{uu}r\equiv0}
% $$
%
% $$
% \kappa_x := \partial_xh+\frac{h}{\rho}\partial_x\rho, \qquad \kappa_y := \partial_yh+\frac{h}{\rho}\partial_y\rho
% $$
%
% $$
% \delta_ur[\phi_i] = \kappa_x\phi_i+h\,\partial_x\phi_i, \qquad \delta_vr[\phi_i] = \kappa_y\phi_i+h\,\partial_y\phi_i
% $$
%
% $$
% \delta^2_{uu}J_{\dot h}[\phi_i,\phi_k] = \frac{1}{\mathcal A}\int\epsilon_{\dot h}^{-2}\big(\kappa_x\phi_i+h\partial_x\phi_i\big)\big(\kappa_x\phi_k+h\partial_x\phi_k\big)\,dx\,dy
% $$
%
% $$
% \delta^2_{vv}J_{\dot h}[\phi_i,\phi_k] = \frac{1}{\mathcal A}\int\epsilon_{\dot h}^{-2}\big(\kappa_y\phi_i+h\partial_y\phi_i\big)\big(\kappa_y\phi_k+h\partial_y\phi_k\big)\,dx\,dy
% $$
%
% $$
% \delta^2_{uv}J_{\dot h}[\phi_i,\phi_k] = \frac{1}{\mathcal A}\int\epsilon_{\dot h}^{-2}\big(\kappa_x\phi_i+h\partial_x\phi_i\big)\big(\kappa_y\phi_k+h\partial_y\phi_k\big)\,dx\,dy
% $$
%
% $$
% J^{qq}_{\dot h} = \frac{1}{\mathcal A}\int \epsilon_{\dot h}^{-2}\;\big(\delta_qr\big)\otimes\big(\delta_qr\big)\;dx\,dy
% $$
%
% $$
% r = \dot h_\mathrm{obs} - a + \left(\partial_xh+\frac{h}{\rho}\partial_x\rho\right)u + h\,\partial_xu + \left(\partial_yh+\frac{h}{\rho}\partial_y\rho\right)v + h\,\partial_yv
% $$
%%




Area=MUA.Area;

duIdu=0;
duIhdot=0;
dvIdv=0;
dvIhdot=0;



if contains(CtrlVar.Inverse.Measurements,"-uv-")
    us=F.ub+F.ud ;
    vs=F.vb+F.vd ;

    uErr=sqrt(spdiags(Meas.usCov));
    usres=(us-Meas.us)./uErr;
    vErr=sqrt(spdiags(Meas.vsCov));
    vsres=(vs-Meas.vs)./vErr;
    duIdu=(MUA.M*usres)./uErr/Area;
    dvIdv=(MUA.M*vsres)./vErr/Area;
    Iuv=full(usres'*MUA.M*usres+vsres'*MUA.M*vsres)/2/Area;

end

% Derivatives of the Ihdot misfit term. This is a bit more complicated as I need to take the derivative of the
% mass-conservation equation with respect to v.
%
% BTW, as the regularization term (R), does not depend on v these I derivatives with respect to v are the derivatives of J
% with respect to v.
%
% There is an interesting added aspect to this, because the mass-conservation involves both v and B, i.e. both the output of
% the forward model (v), and the variable which we are inverting for (B), the 'output' of the inverse model. For this reason,
% we have here a contribution from the misfit term (I) to the quantity <\partial_B J, \phi > which otherwise would only arise
% from the regularization term (R).  Therefore, we must here calculate the derivative of the cost function term with respect
% to both u and B. This addition is only needed when inverting for B while also including the hdot cost function term.
%
if contains(CtrlVar.Inverse.Measurements,"-dhdt-")
    % when including dhdt meas, there is an extra contribution to the RHS of the adjoint system related to
    %
    % $$ \delta_u J_{\dot{h}} $$ and  $$ \delta_v J_{\dot{h}} $$ and 
    %
    %
    [Ihdot,duIhdot,dvIhdot,~]=EvaluateJhdotAndDerivatives([],CtrlVar,MUA,F,BCs,Meas);

end


dJduv=[duIdu(:)+duIhdot(:);dvIdv(:)+dvIhdot(:)];



end





