function varargout=ExplicitEstimationUsingBackwardDifferences(dt,dtm1,dtm2,varargin)
%%
% varargout=ExplicitEstimationUsingBackwardDifferences(dt,dtm1,dtm2,X1,D1,Dm1_1,X2,D2,Dm1_2,...)
%
% Second-order explicit estimate (extrapolation) of fields X at time t(n+1)=t(n)+dt, using backward differences over the two
% previous time steps. (9 Oct 2026)
%
% The inputs are given in triples, one triple for each field:
%
%   X    : X(t(n))
%   D    : (X(t(n))-X(t(n-1)))/dtm1      backward difference over the previous time step,     eg F0.dhdt
%   Dm1  : (X(t(n-1))-X(t(n-2)))/dtm2    backward difference over the time step before that, eg Fm1.dhdt
%
% and
%   dt   : t(n+1)-t(n)       the new time step
%   dtm1 : t(n)-t(n-1)       the time step over which D was calculated   (F0.dtRates)
%   dtm2 : t(n-1)-t(n-2)     the time step over which Dm1 was calculated (Fm1.dtRates)
%
% The estimate is
%
%   X(t(n+1)) = X + dt * ( D + (D-Dm1) * (dtm1+dt)/(dtm2+dtm1) )
%
% ie the mean slope over the new time step is extrapolated linearly from the slopes D and Dm1, which are the time derivatives at the
% midpoints of the two previous time steps. This is equivalent to quadratic extrapolation through X(t(n-2)), X(t(n-1)) and X(t(n)).
% It is exact if X is a quadratic function of time, and it is second-order accurate (local error O(dt^3)), also for variable time steps.
%
% Note: ExplicitEstimation.m implements the variable time step two-step Adams-Bashforth method (AB2), which is correct if its inputs
% are the time derivatives AT t(n) and t(n-1). If the backward differences D and Dm1 are used in AB2 instead, the estimate is only
% first-order accurate; for a constant time step it then reduces to X+dt*dX/dt(t(n)). This function is therefore used for the
% backward differences calculated in UpdateFtimeDerivatives.m. See also TestABExtrapolation.m.
%
% Fallback, node by node:
%
%   - where D and Dm1 are finite, and dtm1 and dtm2 are finite and positive : the second-order estimate above
%   - where D is finite, but not Dm1 (or dtm2 is not finite and positive)   : linear extrapolation, X+dt*D
%   - where D is not finite (or dtm1 is not finite and positive)            : no extrapolation, X
%
% If the number of elements of D or Dm1 does not agree with that of X, these are treated as not available. An empty X is returned as empty.
%
%%

nInputs=nargin-3;
nOutputs=nargout;
varargout=cell(1,nOutputs);

if nInputs~=3*nOutputs
    error('ExplicitEstimationUsingBackwardDifferences:WrongNumberOfInputs',' wrong number of inputs ')
end

isPositiveScalar=@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0 ;
Linear=isPositiveScalar(dtm1) ;
Quadratic=Linear && isPositiveScalar(dtm2) ;

for I=1:nOutputs

    X=varargin{1+3*(I-1)};
    D=varargin{2+3*(I-1)};
    Dm1=varargin{3+3*(I-1)};

    Xest=X;

    if ~isempty(X) && Linear && numel(D)==numel(X)

        iD=isfinite(D);
        Xest(iD)=X(iD)+dt*D(iD);

        if Quadratic && numel(Dm1)==numel(X)
            iDm1=iD & isfinite(Dm1);
            Xest(iDm1)=X(iDm1)+dt*(D(iDm1)+(D(iDm1)-Dm1(iDm1))*(dtm1+dt)/(dtm2+dtm1));
        end

    end

    varargout{I}=Xest;

end

end
