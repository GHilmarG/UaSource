




function [d2taubxdhdh,d2taubydhdh,d2taubxdhdu,d2taubxdhdv,d2taubydhdu,d2taubydhdv]=...
    WeertmanSecondOrderDerivatives(CtrlVar,Dh,ub,vb,C,m)

%% Second-order derivatives of the basal traction for the Weertman sliding law
%
% Returns
%
% $$ \frac{\partial^2 t_{bx}}{\partial h^2} , \quad \frac{\partial^2 t_{by}}{\partial h^2} , \quad
%    \frac{\partial^2 t_{bx}}{\partial h \partial u} , \quad \frac{\partial^2 t_{bx}}{\partial h \partial v} , \quad
%    \frac{\partial^2 t_{by}}{\partial h \partial u} , \quad \frac{\partial^2 t_{by}}{\partial h \partial v} $$
%
% evaluated at the integration points. All outputs beyond the first two are nargout-guarded.
%
%% Why this is a separate function
%
% BasalDrag.m returns first derivatives only, the last two outputs being dtaubxdh and dtaubydh. The second-order
% derivatives are needed for the direct-adjoint Hessian, specifically by
%
%   FBB.m    which needs the h-h derivatives
%   FBuv.m   which needs the mixed h-velocity derivatives
%
% and are supplied here instead.
%
% This is a temporary arrangement. BasalDrag.m is the natural home for these quantities, and when that function is
% rewritten this helper should be folded into it and deleted.
%
%% Weertman only
%
% For the Weertman law
%
% $$ t_{bx} = \mathcal{G}(\Delta h) \, \beta^2 \, u_b , \qquad
%    \beta^2 = (C+C_0)^{-1/m} \, U^{1/m-1} , \qquad U = \sqrt{u_b^2+v_b^2+u_0^2} $$
%
% and $\beta^2$ does not depend on the thickness: the entire thickness dependence sits in the grounding mask. This
% makes every second derivative below a first derivative of BasalDrag.m with one factor replaced:
%
% $$ \frac{\partial^2 t_{bx}}{\partial h^2}
%      = \tilde{\delta}'(\Delta h) \, \beta^2 u_b
%      \qquad \textrm{(i.e. } \partial t_{bx}/\partial h \textrm{ with } \tilde{\delta} \to \tilde{\delta}' ) $$
%
% $$ \frac{\partial^2 t_{bx}}{\partial h \partial u} = \tilde{\delta}(\Delta h) \left ( \beta^2 + D\beta^2 u_b^2 \right ) ,
%    \qquad
%    \frac{\partial^2 t_{bx}}{\partial h \partial v} = \tilde{\delta}(\Delta h) \, D\beta^2 \, u_b v_b $$
%
% which are $\partial t_{bx}/\partial u$ and $\partial t_{bx}/\partial v$ with $\mathcal{G} \to \tilde{\delta}$.
% Compare the Weertman sub-function of BasalDrag.m, lines 493 to 500.
%
% For the effective-pressure dependent laws (Budd, Tsai, Cornford, Umbi) $\beta^2$ depends on $N$, and hence on
% $\Delta h$, so there are further terms and none of the above applies. Those laws are rejected here.
%
%% A note on notation
%
% $\tilde{\delta}$ is the derivative of the SMOOTHED Heaviside function, and is a bounded, ordinary function, not a
% Dirac delta. Its derivative $\tilde{\delta}'$ is therefore an ordinary classical derivative. From
% $\tilde{\delta}(\xi)=2k\mathcal{H}(\xi)(1-\mathcal{H}(\xi))$ it follows that
%
% $$ \tilde{\delta}'(\xi) = 2 k \, \tilde{\delta}(\xi) \, (1-2\mathcal{H}(\xi)) $$
%
% so no further special function is needed.
%
%% Inputs
%
%   Dh        h-h_f , the flotation measure, at the integration points
%   ub,vb     basal velocities at the integration points
%   C,m       slipperiness and sliding exponent at the integration points. C should already have been passed through
%             SmoothFloor.m by the caller, as is done everywhere else.
%
%  see also: BasalDrag.m, FBB.m, FBuv.m, dGeometrydB.m
%
%%

narginchk(6,6)
nargoutchk(2,6)

if ~ismember(string(CtrlVar.SlidingLaw),["W","Weertman"])
    error("WeertmanSecondOrderDerivatives:SlidingLawNotImplemented",...
        ["Second-order derivatives of the basal traction are currently implemented for the Weertman sliding\n" ...
         "law only, but CtrlVar.SlidingLaw=""%s"".\n" ...
         "For the effective-pressure dependent laws beta^2 depends on N, and hence on h, giving further terms."],...
        CtrlVar.SlidingLaw)
end

if CtrlVar.Hh0~=0
    error("WeertmanSecondOrderDerivatives:NonZeroHh0",...
        "A symmetrical Heaviside function is assumed, i.e. CtrlVar.Hh0=0, but CtrlVar.Hh0=%g \n",CtrlVar.Hh0)
end

C0 = CtrlVar.Czero ;
u0 = CtrlVar.SpeedZero ;

U  = sqrt(ub.*ub+vb.*vb+u0*u0) ;
t  = (U./(C+C0)).^(1./m) ;
beta2  = t./U ;                        % as in BasalDrag.m
Dbeta2 = (1./m-1).*t./U.^3 ;           % as in BasalDrag.m , zero for m=1

G      = HeavisideApprox(CtrlVar.kH,Dh,CtrlVar.Hh0) ;
Dtilde = DiracDelta(CtrlVar.kH,Dh,CtrlVar.Hh0) ;        % regularised delta, not a Dirac delta
Dprime = 2*CtrlVar.kH.*Dtilde.*(1-2*G) ;                % its ordinary derivative

%% second derivatives with respect to the thickness

d2taubxdhdh = Dprime.*beta2.*ub ;
d2taubydhdh = Dprime.*beta2.*vb ;

%% mixed thickness-velocity derivatives
%
% These are the velocity derivatives of BasalDrag.m with the grounding mask replaced by the regularised delta.

if nargout>2 ; d2taubxdhdu = Dtilde.*( beta2 + Dbeta2.*ub.*ub ) ; end
if nargout>3 ; d2taubxdhdv = Dtilde.*Dbeta2.*ub.*vb ;             end
if nargout>4 ; d2taubydhdu = Dtilde.*Dbeta2.*vb.*ub ;             end
if nargout>5 ; d2taubydhdv = Dtilde.*( beta2 + Dbeta2.*vb.*vb ) ; end

end
