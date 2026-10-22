function out=BoxTransform(Mode,x,Box)

%%
%
%   v = BoxTransform("forward" ,p,Box)     maps  p in (lb,ub)  ->  v in R
%   p = BoxTransform("inverse" ,v,Box)     maps  v in R        ->  p in (lb,ub)
%   D = BoxTransform("jacobian",v,Box)     returns dp/dv evaluated at v
%
% Change of variables used to eliminate two-sided box constraints. Box is built once by BoxTransformSetup.m, which
% documents the two available maps and the reason the softplus box is the default.
%
%
%% The softplus box
%
%       g(x,s) = s*log(1+exp(x/s))
%       p      = v - g(v-ub,s) + g(lb-v,s)
%       dp/dv  = 1 - sigma((v-ub)/s) - sigma((lb-v)/s)
%
% dp/dv lies in (0,1): it tends to one in the interior, so the transformed problem is essentially the original one
% there, and tends to zero at either bound, so the bounds hold for every v. Crucially it can never exceed one, so a
% gradient can never be amplified by the transformation.
%
% The forward map has no closed form and is obtained by Newton iteration. Since p is close to v this converges in a
% couple of steps for interior variables, and it is only ever needed once, when F2p maps the starting fields.
%
%
%% The anchored logistic
%
%       p      = lb + (ub-lb)*sigma( logit(r0) + (v-p0)/D0 ) ,  D0 = dLo0*dHi0/W
%       dp/dv  = dLo*dHi / (dLo0*dHi0)
%
% Exactly one at the anchor, but unbounded above thereafter. Kept as an option.
%
%
%% A note on where the Jacobian belongs
%
% dJ/dv = (dp/dv) * dJ/dp is a relation between l2 gradients. Where the transformation is combined with a Riesz map, the
% chain rule must therefore be applied to the l2 gradient BEFORE the Riesz map. Applying it afterwards would require
% G\(D.*(G*dJdp)) instead of D.*dJdp, i.e. an extra matrix-vector product and an extra solve with G on every
% cost-function evaluation. See JGH.m .
%
%%

narginchk(3,3)

x=x(:) ;

if ~isstruct(Box) || ~isfield(Box,'Type')
    error('BoxTransform:BadBox','Box must be the structure returned by BoxTransformSetup.')
end

if numel(x)~=numel(Box.Ok)
    error('BoxTransform:SizeMismatch','x has %i elements but the transformation was set up for %i.',...
        numel(x),numel(Box.Ok))
end

Ok=Box.Ok ;
out=x ;

switch Box.Type

    case "softplus"

        lb=Box.Lower(Ok) ; ub=Box.Upper(Ok) ;
        s=Box.Width_s(Ok) ;

        switch lower(string(Mode))

            case "forward"

                dF=Box.DistanceFloor(Ok) ;
                q=min(max(x(Ok),lb+dF),ub-dF) ;
                out(Ok)=SoftBoxForward(q,lb,ub,s) ;

            case "inverse"

                out(Ok)=SoftBoxInverse(x(Ok),lb,ub,s) ;

            case "jacobian"

                out=ones(size(x)) ;
                out(Ok)=SoftBoxJacobian(x(Ok),lb,ub,s) ;

            otherwise

                error('BoxTransform:UnknownMode','Mode must be "forward", "inverse" or "jacobian", but is "%s".',string(Mode))

        end

    case "logistic"

        switch lower(string(Mode))

            case "forward"

                dF=Box.DistanceFloor(Ok) ;
                q=min(max(x(Ok),Box.Lower(Ok)+dF),Box.Upper(Ok)-dF) ;
                dLo=q-Box.Lower(Ok) ;
                dHi=Box.Upper(Ok)-q ;
                out(Ok)=Box.pAnchor(Ok)+Box.Scale(Ok).*(log(dLo)-log(dHi)-Box.LogitAnchor(Ok)) ;

            case "inverse"

                [dLo,dHi]=LogisticDistances(x,Box,Ok) ;
                p=Box.Lower(Ok)+dLo ;
                Near=dHi<dLo ;
                pHi=Box.Upper(Ok)-dHi ;
                p(Near)=pHi(Near) ;
                out(Ok)=p ;

            case "jacobian"

                [dLo,dHi]=LogisticDistances(x,Box,Ok) ;
                out=ones(size(x)) ;
                out(Ok)=dLo.*dHi./Box.Width(Ok)./Box.Scale(Ok) ;

            otherwise

                error('BoxTransform:UnknownMode','Mode must be "forward", "inverse" or "jacobian", but is "%s".',string(Mode))

        end

    otherwise

        error('BoxTransform:UnknownType','Box.Type must be "softplus" or "logistic", but is "%s".',Box.Type)

end

end

%%%%%%%%%%%%%%%% softplus box %%%%%%%%%%%%%%%%

function p=SoftBoxForward(q,lb,ub,s)

% p -> v by Newton iteration on  SoftBoxInverse(v)=q .
%
% dp/dv is in (0,1) and p(v) is strictly increasing, so the iteration is well behaved. Since p is close to v away from
% the bounds, starting from v=q converges in a couple of steps for most variables.

p=q ;   % here p is being used as v, the iterate

Tol=1e-12*max(abs(q),1) ;
Residual=inf(size(q)) ;

for Iteration=1:200

    Residual=SoftBoxInverse(p,lb,ub,s)-q ;

    if all(abs(Residual)<Tol) ; break ; end

    D=SoftBoxJacobian(p,lb,ub,s) ;
    D=max(D,1e-14) ;

    Step=-Residual./D ;

    % cap the step at a generous multiple of the local barrier width, so that a variable starting on a bound climbs out
    % steadily rather than overshooting to where the exponentials underflow
    Cap=50*s ;
    Step=max(min(Step,Cap),-Cap) ;

    p=p+Step ;

end

if any(abs(Residual)>=Tol)
    warning('BoxTransform:ForwardNotConverged',...
        'The Newton iteration for the forward map did not converge for %i of %i variables, max residual %g.',...
        sum(abs(Residual)>=Tol),numel(q),max(abs(Residual)))
end

end

%%

function p=SoftBoxInverse(v,lb,ub,s)

% p = v - g(v-ub,s) + g(lb-v,s) , with g the softplus.
%
% Three algebraically identical forms are used, one per regime, to avoid cancellation. In the interior p is close to v,
% so v plus two small corrections is accurate. Above ub, or below lb, that form would subtract two large nearly equal
% numbers, and the identity  v - g(v-ub,s) = ub - g(ub-v,s)  is used instead.

p=zeros(size(v)) ;

iHi = v>=ub ;
iLo = v<=lb ;
iMid= ~iHi & ~iLo ;

p(iMid)=v(iMid)-SoftPlus(v(iMid)-ub(iMid),s(iMid))+SoftPlus(lb(iMid)-v(iMid),s(iMid)) ;
p(iHi) =ub(iHi)-SoftPlus(ub(iHi)-v(iHi),s(iHi))  +SoftPlus(lb(iHi)-v(iHi),s(iHi)) ;
p(iLo) =lb(iLo)+SoftPlus(v(iLo)-lb(iLo),s(iLo))  -SoftPlus(v(iLo)-ub(iLo),s(iLo)) ;

% The limits are approached but never crossed analytically; this only removes the last rounding error.
p=min(max(p,lb),ub) ;

end

%%

function D=SoftBoxJacobian(v,lb,ub,s)

D=1-Sigma((v-ub)./s)-Sigma((lb-v)./s) ;

end

%%

function y=SoftPlus(x,s)

% s*log(1+exp(x/s)), evaluated so that it does not overflow for large x/s

z=x./s ;
y=zeros(size(z)) ;

i=z>0 ;
y(i)=x(i)+s(i).*log1p(exp(-z(i))) ;

i=~i ;
y(i)=s(i).*log1p(exp(z(i))) ;

end

%%

function y=Sigma(z)

% the logistic function, written with tanh so that it does not overflow

y=0.5*(1+tanh(z/2)) ;

end

%%%%%%%%%%%%%%%% anchored logistic %%%%%%%%%%%%%%%%

function [dLo,dHi]=LogisticDistances(v,Box,Ok)

% Distances to the two bounds at v, i.e. dLo=p-lb and dHi=ub-p, with dLo+dHi=W. The two branches avoid overflow and
% keep full precision in whichever distance is the small one.

z=Box.LogitAnchor(Ok)+(v(Ok)-Box.pAnchor(Ok))./Box.Scale(Ok) ;
W=Box.Width(Ok) ;

dLo=zeros(size(z)) ;
dHi=zeros(size(z)) ;

i=z>=0 ;
dHi(i)=W(i)./(1+exp(z(i))) ;
dLo(i)=W(i)-dHi(i) ;

i=~i ;
dLo(i)=W(i)./(1+exp(-z(i))) ;
dHi(i)=W(i)-dLo(i) ;

end
