function [h,dPhi,d2Phi,x,Active]=ThicknessFloor(CtrlVar,v,Mode)

%% Smooth floor on the ice thickness
%
%   [h,dPhi,d2Phi,x,Active]=ThicknessFloor(CtrlVar,x,"forward")
%   [h,dPhi,d2Phi,x,Active]=ThicknessFloor(CtrlVar,h,"inverse")
%
% Maps a thickness $x$ to a thickness $h$ that is never smaller than the minimum thickness $h_{\min}$=CtrlVar.ThickMin. The map is
%
% $$ h = \Phi(x) = h_{\min} + w \, \ln \left( 1 + \exp \left( \frac{x-h_{\min}}{w} \right) \right) $$
%
% where the width, $w$, is
%
% $$ w = \mathrm{CtrlVar.ThickMinWidthRelative} \times h_{\min} $$
%
% This is the smooth floor, or softplus, function. It is used in Calc_bh_From_sBS.m, where $x=h_c$ is the thickness solving the geometrical
% closure and $h=\Phi(h_c)$ the thickness returned, and in dGeometrydB.m, where the derivatives of the closure are modified accordingly.
%
%% Properties
%
% * $\Phi(x) > h_{\min}$ for all $x$, and $\Phi(x) \rightarrow h_{\min}$ as $x \rightarrow -\infty$.
%
% * $\Phi(x) \approx x$ for $x \gg h_{\min}$. The difference is $w \exp(-(x-h_{\min})/w)$, which is below $10^{-13} \, w$ for $x>h_{\min}+30 w$. Where
%   $x > h_{\min}+40 w$ the thickness is returned exactly unchanged, so that thick ice is not affected at all.
%
% * $\Phi(h_{\min}) = h_{\min} + w \ln 2$, i.e. about $0.035$ m for $h_{\min}=1$ m and the default width.
%
% * The limit $w \rightarrow 0$ is $\max(x,h_{\min})$. A hard maximum is not used because it has a kink, and is not invertible.
%
% * The first and second derivatives are
%
% $$ \Phi'(x) = \frac{1}{1+\exp(-(x-h_{\min})/w)} , \qquad \Phi''(x) = \frac{\Phi'(x) \, (1-\Phi'(x))}{w} $$
%
% * The inverse, which exists for $h>h_{\min}$, is
%
% $$ x = \Phi^{-1}(h) = h_{\min} + w \ln \left( \exp \left( \frac{h-h_{\min}}{w} \right) - 1 \right) $$
%
%% Inputs and outputs
%
% Mode="forward": the input, v, is the thickness $x$ and the output h is $\Phi(x)$.
%
% Mode="inverse": the input, v, is the (floored) thickness $h$ and the unfloored thickness $x=\Phi^{-1}(h)$ is returned as the fourth output. Values
% $h \le h_{\min}$ can not be the result of the floor, and are treated as $h$ one machine epsilon above $h_{\min}$.
%
% Because of round-off, $\Phi(x)$ is equal to $h_{\min}$ to machine precision for $x < h_{\min}-35 w$, and $x$ can then not be recovered from $h$. The inverse then
% returns a value close to $h_{\min}-36 w$. This has no consequence for the derivatives, which are zero to machine precision for such $x$ in either case.
%
% The derivatives dPhi and d2Phi are always evaluated at $x$, i.e. at the unfloored thickness. Outputs:
%
%   h        the floored thickness, $\Phi(x)$ (in "inverse" mode this is the input).
%   dPhi     $\Phi'(x)$
%   d2Phi    $\Phi''(x)$
%   x        the unfloored thickness (in "forward" mode this is the input).
%   Active   false if the floor is switched off, in which case h=x, dPhi=1 and d2Phi=0.
%
%% Options
%
%   CtrlVar.ThickMin               the minimum thickness, $h_{\min}$.
%   CtrlVar.ThickMinWidthRelative  the width of the floor in units of ThickMin (default 0.05). If zero, or if ThickMin is zero, there is no floor.
%
% If CtrlVar does not have the field ThickMinWidthRelative, for example a CtrlVar read from an old restart file, the default value 0.05 is
% used. This is the same value as is set in Ua2D_DefaultParameters.m.
%
%  see also: Calc_bh_From_sBS.m, dGeometrydB.m
%
%%

narginchk(3,3)
nargoutchk(1,5)

hmin=CtrlVar.ThickMin ;

if isfield(CtrlVar,"ThickMinWidthRelative") && ~isempty(CtrlVar.ThickMinWidthRelative)
    WidthRelative=CtrlVar.ThickMinWidthRelative ;
else
    WidthRelative=0.05 ;
end

w=WidthRelative*hmin ;

Active=(hmin>0) && (w>0) ;

if ~Active
    h=v ; x=v ;
    dPhi=ones(size(v)) ;
    d2Phi=zeros(size(v)) ;
    return
end

zThick=40 ;     % beyond this many widths above the minimum thickness the floor is exactly the identity, to machine precision

switch lower(string(Mode))

    case "forward"

        x=v ;
        z=(x-hmin)/w ;
        h=hmin+w*(max(z,0)+log1p(exp(-abs(z)))) ;
        Thick=z>zThick ;
        h(Thick)=x(Thick) ;

    case "inverse"

        h=v ;
        y=max((h-hmin)/w,eps) ;             % y = softplus(z)
        z=y+log(-expm1(-y)) ;               % z = inverse of softplus
        x=hmin+w*z ;
        Thick=y>zThick ;
        x(Thick)=h(Thick) ;
        z(Thick)=(x(Thick)-hmin)/w ;

    otherwise

        error("ThicknessFloor:Mode","Mode must be ""forward"" or ""inverse"".")

end

% derivatives at x, written in a form that does not overflow
E=exp(-abs(z)) ;
dPhi=1./(1+E) ;
Negative=z<0 ;
dPhi(Negative)=E(Negative)./(1+E(Negative)) ;
d2Phi=E./(1+E).^2/w ;

end
