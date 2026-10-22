function Box=BoxTransformSetup(p,plb,pub,Type,WidthFraction)

%%
%
%   Box=BoxTransformSetup(p,plb,pub)
%   Box=BoxTransformSetup(p,plb,pub,Type,WidthFraction)
%
% Builds a change of variables that eliminates two-sided box constraints. The resulting structure is passed to
% BoxTransform.m and is normally stored on MUA so that the same transformation is used everywhere and for the whole run.
%
%   Type            "softplus" (default) or "logistic"
%   WidthFraction   shape parameter, see below. Default 0.1 for softplus, unused for logistic.
%
%
%% Why "softplus" is the default
%
% The obvious choice, p = lb + (ub-lb)*sigma(u), is a poor one whenever a variable sits close to a bound but belongs far
% from it. Its Jacobian is
%
%       dp/du = dLo*dHi / (dLo0*dHi0)
%
% which is UNBOUNDED ABOVE: a variable starting 0.1 m from its bound and moving to 100 m from it has its gradient
% amplified by a factor of a thousand on the way out. On a Hofsjokull bed inversion that produced decrements jumping by
% seven orders of magnitude between iterations, CG correction ratios of order 200, and repeated zero steps. Anchoring
% fixes the scaling at the starting point but cannot fix the shape of the function.
%
% The softplus box instead behaves like the identity in the interior and bends only close to a bound:
%
%       g(x,s) = s*log(1+exp(x/s))                       softplus of width s
%       p      = v - g(v-ub,sHi) + g(lb-v,sLo)
%       dp/dv  = 1 - sigma((v-ub)/sHi) - sigma((lb-v)/sLo)
%
% The Jacobian lies in (0,1). It CANNOT amplify, it tends to one in the interior, and it tends to zero at either bound,
% so the bounds are still enforced for every v. On the same fields the Jacobian ranges over [0.5,1] with a median of
% exactly 1, against [1,2000] for the logistic.
%
% Because p is close to v, the new variable never runs away either: reaching within delta of a bound needs only
% v-ub = s*log(s/delta), so about 2 m for s=0.5 m and delta=1 cm.
%
%
%% The width
%
% A single width per variable is used for both barriers,
%
%       s = WidthFraction * min( p0-lb , ub-p0 )
%
% i.e. it is set by the NEARER of the two bounds. At the starting point the near barrier then contributes
% sigma(-1/WidthFraction) to the Jacobian, about 4.5e-5 for WidthFraction=0.1, and the far one contributes less, so
%
%       dp/dv = 1 - sigma(-1/WidthFraction) - (something smaller)
%
% which is the same for every variable to within a few parts in 1e5. The Jacobian is therefore uniform at the start by
% construction, as the anchored logistic also achieved, but now with the guarantee that it can never exceed one.
%
% Taking the minimum is what keeps a non-binding bound harmless, and it has to be done this way. Scaling each barrier by
% its own distance instead, sLo = WidthFraction*(p0-lb), fails badly for a distant bound: with lb 1e9 m below the bed,
% sLo would be 1e8 m, and although that barrier is flat there, its residual value s*log(1+exp(-1/WidthFraction)) is then
% about 4500 m, a large constant offset between v and p. With s set by the nearer bound, the far barrier's contribution
% is exp(-distance/s) with distance/s of order 1e7, i.e. exactly zero, and the transformation is identical whether the
% flotation bound is 9 km below the bed or 1e12 m below it.
%
%
%% Untransformed components
%
% Components whose box is degenerate (ub=lb, as produced by ubB=max(lbB,ubB) in F2p where the flotation bound exceeds
% the surface bound) or unbounded (lb=-inf or ub=inf) are left as they are: the map is the identity there and the
% Jacobian is one. Box.Ok records which components are actually transformed.
%
%%

narginchk(3,5)

if nargin<4 || isempty(Type)          ; Type="softplus" ; end
if nargin<5 || isempty(WidthFraction) ; WidthFraction=0.1 ; end

Type=lower(string(Type)) ;

p=p(:) ; plb=plb(:) ; pub=pub(:) ;

if isempty(plb) || isempty(pub)
    error('BoxTransformSetup:NoBounds','plb and pub must not be empty.')
end

W=pub-plb ;

Box.Type=Type ;
Box.Ok = isfinite(W) & W>0 & isfinite(p) ;
Box.Lower=plb ;
Box.Upper=pub ;
Box.Width=W ;

% Floor at rounding level, only so that the anchor distances are strictly positive when p sits exactly on a bound, which
% is the normal situation when restarting from a run that used the bounded optimiser. This is deliberately NOT a
% fraction of the box width: with a distant flotation bound W can be 1e9 m, and a floor of 1e-6*W would be 1000 m,
% clamping perfectly ordinary values.
Box.DistanceFloor=4*eps(W) ;

pAnchor=p ;
pAnchor(Box.Ok)=min(max(p(Box.Ok),plb(Box.Ok)+Box.DistanceFloor(Box.Ok)),pub(Box.Ok)-Box.DistanceFloor(Box.Ok)) ;
Box.pAnchor=pAnchor ;

dLo0=zeros(size(p)) ; dHi0=zeros(size(p)) ;
dLo0(Box.Ok)=pAnchor(Box.Ok)-plb(Box.Ok) ;
dHi0(Box.Ok)=pub(Box.Ok)-pAnchor(Box.Ok) ;

Box.rAnchor=ones(size(p)) ;
Box.rAnchor(Box.Ok)=dLo0(Box.Ok)./W(Box.Ok) ;

switch Type

    case "softplus"

        Box.WidthFraction=WidthFraction ;

        s=ones(size(p)) ;
        s(Box.Ok)=WidthFraction*min(dLo0(Box.Ok),dHi0(Box.Ok)) ;

        % A variable sitting exactly on a bound would otherwise get a barrier of zero width, i.e. an infinitely sharp
        % one. Floor at a small fraction of the typical width so that it behaves like a variable that is merely very
        % close to its bound. Its dp/dv at the anchor is then well below one, which is the intended behaviour.
        Box.WidthFloorFraction=1e-4 ;
        s=FloorAtFractionOfMedian(s,Box.Ok,Box.WidthFloorFraction) ;

        Box.Width_s=s ;

        % The two barriers must leave a positive Jacobian everywhere. With s set by the nearer distance this holds for
        % WidthFraction below about 0.4, but check rather than assume.
        D0=BoxTransform("jacobian",BoxTransform("forward",pAnchor,Box),Box) ;
        if any(D0(Box.Ok)<=0)
            error('BoxTransformSetup:BarriersOverlap',...
                ['The two softplus barriers overlap: dp/dv at the starting point is %g for at least one variable.\n',...
                'Reduce CtrlVar.Inverse.BoxTransformWidthFraction, currently %g.'],min(D0(Box.Ok)),WidthFraction)
        end
        Box.JacobianAtAnchor=D0 ;

    case "logistic"

        % Anchored logistic, kept as an option. See the note above on why it is not the default.
        Box.LogitAnchor=zeros(size(p)) ;
        Box.Scale=ones(size(p)) ;

        Box.LogitAnchor(Box.Ok)=log(dLo0(Box.Ok))-log(dHi0(Box.Ok)) ;

        D0=dLo0.*dHi0./W ;
        Box.ScaleFloorFraction=1e-3 ;
        ScaleFloored=FloorAtFractionOfMedian(D0,Box.Ok,Box.ScaleFloorFraction) ;
        Box.Scale(Box.Ok)=ScaleFloored(Box.Ok) ;
        Box.JacobianAtAnchor=ones(size(p)) ;

    otherwise

        error('BoxTransformSetup:UnknownType','Type must be "softplus" or "logistic", but is "%s".',Type)

end

end

%%

function x=FloorAtFractionOfMedian(x,Ok,Fraction)

% Floor the transformable entries of x at Fraction times their median, leaving the rest at one.

xOk=x(Ok) ;
xTypical=median(xOk(xOk>0)) ;

if isempty(xTypical) || ~isfinite(xTypical) || xTypical<=0 ; xTypical=1 ; end

x(Ok)=max(x(Ok),Fraction*xTypical) ;

end
