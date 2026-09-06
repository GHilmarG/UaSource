

function [d1,cgInfo]=NewConjugatedGradMetric(mdJdp,mdJdplast,d0,G,CtrlVar,cgInfo)

%%
%
%   [d1,RunInfo]=NewConjugatedGradMetric(mdJdp,mdJdplast,d0,G,CtrlVar,RunInfo)
%
% Conjugated-gradient update of the search direction for a general metric.
%
% Metric-aware version of NewConjugatedGrad.m . Every inner product is taken
% with respect to the metric matrix G, i.e.
%
%       <x,y>_G = x' G y        ||x||_G = sqrt( x' G x )
%
% rather than the Euclidean products x'*y and sqrt(x'*x). Setting G=[] (or G=I)
% recovers the Euclidean update.
%
%
%% Sign convention  (differs from NewConjugatedGrad.m, see note at the end)
%
% Write g for the gradient after the Riesz map has been applied (see below).
% The inputs are the steepest-DESCENT directions
%
%       mdJdp     = -dJdp        (current,  "minus dJdp")
%       mdJdplast = -dJdplast    (previous)
%
% and the output is the new search DESCENT direction, so that the line search is
%
%       minimise  J(p0 + gamma*d1)   over   gamma > 0 .
%
% The returned d1 is fed back in as d0 on the next call, unaltered and with no
% change of sign. Do not feed d1 back in as mdJdplast.
%
%
%% Inputs
%
%   mdJdp           current steepest-descent direction, -dJdp
%   mdJdplast       previous steepest-descent direction, -dJdplast
%   d0              previous search direction, i.e. the d1 returned by the
%                   previous call to this routine
%   G               metric matrix, must be symmetric positive definite. Given as
%                       []                : Euclidean metric, G=I
%                       a (sparse) matrix : products formed as G*x
%                       a function handle : G(x) must return G*x
%   CtrlVar         uses  CtrlVar.ConjugatedGradientsUpdate             ('FR'|'PR'|'HS'|'DY')
%                         CtrlVar.ConjugatedGradientsSufficientDescent   (default 0.1)
%                         CtrlVar.Inverse.InfoLevel
%   cgInfo          uses/updates cgInfo.ConjGradUpdate
%                   returns      cgInfo.ConjGradAngle , cgInfo.ddAngle
%
%
%% Restarts
%
% The CG update is discarded and d1 reset to the steepest-descent direction s1
% if any of the following holds:
%
%   1) ||s1||_G or ||s0||_G is at roundoff level
%   2) teta is not finite
%   3) teta<0 . For Polak-Ribiere this is the PR+ variant, max(teta,0), which is
%      not merely a safeguard: plain PR can fail to converge on smooth functions
%      (Powell 1984), whereas PR+ is globally convergent under a suitable line
%      search (Gilbert & Nocedal 1992)
%   4) the sufficient-descent condition <s1,d1>_G >= c ||s1||_G^2 is violated,
%      with c = CtrlVar.ConjugatedGradientsSufficientDescent
%
% ddAngle, the angle between successive steepest-descent directions, is computed
% and returned in cgInfo as a diagnostic but is NOT used to trigger a restart.
% NewConjugatedGrad.m restarted on abs(ddAngle)<0.5 degrees, which fires only
% when s1 and s0 are nearly parallel, not when the angle departs from 90 degrees
% as the comment there claimed, and in practice never fired. The scale-free
% version of that test is Powell's, |<s1,s0>_G| >= nu ||s1||_G^2 with nu=0.1
% (Nocedal & Wright 5.52); note the norm ratio in the denominator, which a plain
% angle test cannot express. It is not implemented here, 3) and 4) being enough
% for the convergence theory.
%
% Output
%
%   d1              new search direction, a descent direction
%
%
%% Gradients versus directions
%
% In a non-Euclidean metric one must distinguish between
%
%   g : the gradient returned by the adjoint calculation. This is a linear
%       functional (a co-vector), defined through  DJ(p)[dp] = g' dp .
%
%   s : the corresponding steepest-descent direction (a vector), i.e. the Riesz
%       representer of -g with respect to <.,.>_G :
%
%           <s,v>_G = -g' v   for all v        =>        G s = -g
%
% This routine assumes the Riesz map has already been applied by the caller, so
% that mdJdp is a direction and not a raw l2 co-vector. Each l2 product x'*y of
% the original routine is then replaced by x'*G*y, which is what is done below.
%
% The metric R used for the Riesz map and the metric G used here should be the
% same matrix. If they differ, expanding e.g. Fletcher-Reeves gives
%
%       ||s1||_G^2 / ||s0||_G^2 = g1' inv(R) G inv(R) g1 / ( g0' inv(R) G inv(R) g0 )
%
% so the effective preconditioner becomes the sandwich inv(R) G inv(R). With
% R=G this collapses to the usual preconditioned form g1' inv(G) g1 / g0' inv(G) g0 .
%
%
%% Relation to NewConjugatedGrad.m
%
% NewConjugatedGrad.m negates its input internally (d=-dJdescent) and therefore
% returns a direction of the opposite sign to its input. This routine does not
% negate: descent goes in, descent comes out. If you are switching over from
% NewConjugatedGrad.m, the step in the line search changes from
%
%       p - gamma*gradJ1        to        p + gamma*d1 .
%
% teta is unchanged by the switch for FR and PR, these being ratios of terms
% quadratic in the inputs. For HS and DY it is not: the denominator in
% NewConjugatedGrad.m carries the wrong sign, see the note further down and
% TestNewConjugatedGradMetric.m . That has been corrected here.
%
%%

persistent iCount

if nargin<4 ; G=[] ; end
if nargin<5 ; CtrlVar=[] ; end

if isempty(iCount) ; iCount=0 ; end

iCount=iCount+1;

if nargin<6 || ~isstruct(cgInfo) ; cgInfo=struct ; end

if ~isfield(cgInfo,"NumberOfConjGradUpdatesWithoutReset")
    cgInfo.NumberOfConjGradUpdatesWithoutReset=0;
end

cgInfo.NumberOfConjGradUpdatesWithoutReset=cgInfo.NumberOfConjGradUpdatesWithoutReset+1; 

if ~isfield(cgInfo,"ConjGradAngle")
    cgInfo.ConjGradAngle=nan;
end

if ~isfield(cgInfo,"ddAngle")
    cgInfo.ddAngle=nan;
end


% Sufficient-descent parameter c. c=0 reduces the test to the bare descent
% requirement <s1,d1>_G>0 ; c=1 forces a restart at every iteration.
if isempty(CtrlVar) || ~isfield(CtrlVar,'ConjugatedGradientsSufficientDescent') ...
        || isempty(CtrlVar.ConjugatedGradientsSufficientDescent)
    CtrlVar.ConjugatedGradientsSufficientDescent=0.1 ;
end

if ~isfield(CtrlVar,'ConjugatedGradientsUpdate')
    CtrlVar.ConjugatedGradientsUpdate='FR';
end

if ~isfield(CtrlVar,'Inverse') || ~isfield(CtrlVar.Inverse,'InfoLevel')
    CtrlVar.Inverse.InfoLevel=0;
end

 

%% define the metric operator  x -> G x

if isempty(G)
    Gmul=@(x) x ;
elseif isa(G,'function_handle')
    Gmul=G ;
else
    Gmul=@(x) G*x ;
end

%% steepest-descent directions and their G-products
%
% s1 and s0 are the current and previous steepest-descent directions. Only three
% matrix-vector products with G are formed (Gs0, Gs1, Gd1), everything is reused.

s1=mdJdp(:) ;         % = -dJdp      , current  steepest descent
s0=mdJdplast(:) ;     % = -dJdplast  , previous steepest descent
d0=d0(:) ;            % previous search direction

Gs0=Gmul(s0) ;
Gs1=Gmul(s1) ;

n0Sqr=max(s0'*Gs0,0) ;   % ||s0||_G^2 , max(.,0) guards against roundoff
n1Sqr=max(s1'*Gs1,0) ;   % ||s1||_G^2
n0=sqrt(n0Sqr) ;
n1=sqrt(n1Sqr) ;

s1Gs0=s1'*Gs0 ;          % <s1,s0>_G

% angle between the current and the previous steepest-descent directions,
% measured in the G metric
if n0>0 && n1>0
    dd=s1Gs0/(n1*n0) ;
    dd=min(max(dd,-1),1) ;   % protect acosd against roundoff
else
    dd=0 ;
end
ddAngle=acosd(dd);  % diagnostic only, see the note on restarts in the header

%%

if n1<eps || n0<eps
    
    % return the (G-metric) steepest-descent direction
    ConjGradAngle=0 ; teta=0; d1=s1;
    iCount=0;
    cgInfo.NumberOfConjGradUpdatesWithoutReset=0; 
 
    
else
    
    % seems that Polak-Ribiere is the best update
    %
    % NOTE on HS and DY: the denominator is <s0-s1,d0>_G , which is the G-metric
    % form of the (g1-g0)'p0 denominator of Nocedal & Wright (5.46) and (5.49).
    % NewConjugatedGrad.m had this with the opposite sign, which made teta
    % negative whenever the Polak-Ribiere numerator was positive, so that the
    % teta<0 reset below fired and HS and DY silently degraded to steepest
    % descent. See TestNewConjugatedGradMetric.m .
    %
    % With d0 a descent direction and the Wolfe conditions satisfied we have
    % (g1-g0)'d0 > 0, and since s1-s0 = -(g1-g0) the denominator below is then
    % positive, which is what guarantees teta>0 for Dai-Yuan.
    
    den=d0'*(Gs0-Gs1) ;   % <s0-s1,d0>_G , used by HS and DY
    
    switch upper(CtrlVar.ConjugatedGradientsUpdate)
        
        case 'FR'
            teta=n1Sqr/n0Sqr;                    % Fletcher-Reeves
            
        case 'PR'
            teta=(n1Sqr-s1Gs0)/n0Sqr;            % Polak-Ribiere
            
        case 'HS'
            teta=(n1Sqr-s1Gs0)/den;              % Hestenes-Stiefel
            
        case 'DY'
            teta=n1Sqr/den;                      % Dai-Yuan
            
        otherwise
            error('NewConjugatedGradMetric:CaseNotFound','CtrlVar.ConjugatedGradientsUpdate=%s not recognized',...
                CtrlVar.ConjugatedGradientsUpdate)
    end

    if ~isfinite(teta)

        if CtrlVar.Inverse.InfoLevel>=2
            fprintf(' resetting conjugated gradients because teta is not finite. \n ')
        end
        iCount=0;
        cgInfo.NumberOfConjGradUpdatesWithoutReset=0;

        teta=0;

    end
    
    if teta<0
        
        if CtrlVar.Inverse.InfoLevel>=2
            fprintf(' resetting conjugated gradients because teta=%g<0. \n ',teta)
        end
        iCount=0;
        cgInfo.NumberOfConjGradUpdatesWithoutReset=0; 
   
        teta=0;
        
    end
    
    d1=s1 + teta * d0 ; % new search direction
    
    % Angle between the current steepest-descent and the conj-grad search
    % direction, measured in the G metric.
    %
    % Note that DJ(p)[d1] = g' d1 = -<s1,d1>_G , so d1 is a descent direction
    % if and only if <s1,d1>_G > 0, i.e. if and only if ConjGradAngle < 90 deg.
    
    Gd1=Gmul(d1) ;
    nd1=sqrt(max(d1'*Gd1,0)) ;
    s1Gd1=s1'*Gd1 ;
    
    if n1>0 && nd1>0
        cosAngle=min(max(s1Gd1/(n1*nd1),-1),1) ;
    else
        cosAngle=1 ;
    end
    ConjGradAngle=acosd(cosAngle);
    
    % Sufficient-descent safeguard.
    %
    % In exact arithmetic with an exact line search <s1,d0>_G=0, so that
    % <s1,d1>_G = ||s1||_G^2 exactly. The ratio <s1,d1>_G/||s1||_G^2 is therefore
    % a dimensionless measure of how far the update has degraded, equal to one in
    % the ideal case. Requiring only <s1,d1>_G>0 is weak: a d1 almost G-orthogonal
    % to s1 passes the test while being of little use as a search direction.
    
    if s1Gd1 < CtrlVar.ConjugatedGradientsSufficientDescent*n1Sqr
        if CtrlVar.Inverse.InfoLevel>=2
            fprintf(' resetting conjugated gradients because the sufficient-descent condition is violated: \n')
            fprintf('  <s1,d1>_G/||s1||_G^2=%g < CtrlVar.ConjugatedGradientsSufficientDescent=%g , angle=%g deg. \n',...
                s1Gd1/n1Sqr,CtrlVar.ConjugatedGradientsSufficientDescent,ConjGradAngle)
        end
        iCount=0;
        cgInfo.NumberOfConjGradUpdatesWithoutReset=0; 
      
        teta=0;
        d1=s1;
        ConjGradAngle=0;
    end
    
end

%%

cgInfo.ConjGradAngle=ConjGradAngle;
cgInfo.ddAngle=ddAngle;
cgInfo.teta=teta;

if CtrlVar.Inverse.InfoLevel>=2
    fprintf(' Conj. Grad update %s # %-i \n ',CtrlVar.ConjugatedGradientsUpdate,iCount)
    fprintf(' teta=%g \n',teta)
    fprintf('angle between current and previous steepest descent directions is %-20.10g degrees \n',ddAngle)
    fprintf('     angle between conj. grad. and steepest descent directions is %-20.10g degrees \n',ConjGradAngle)
end

end
